function Results = bst_inverse_dirac(HeadModel, OPTIONS)
% BST_INVERSE_DIRAC: Minimum-norm inverse in the Dirac eigenmode basis (built incrementally).
%
% USAGE:  Results = bst_inverse_dirac(HeadModel, OPTIONS)
%
% DESCRIPTION:
%     Re-implements the whitened minimum-norm method of BST_INVERSE_LINEAR_2018,
%     but with the source side expressed in the curvature-aware Dirac eigenbasis
%     produced by BST_DIRAC. The sensor side is unchanged, so every sensor-space
%     step (noise whitening, SNR regularization) is identical to the canonical MNE;
%     only the forward operator and the final source reconstruction differ.
%
%     This file is built STAGE BY STAGE so each step can be tested against the
%     canonical pipeline before the next is added:
%       STAGE 1 -- Noise-covariance regularization + whitener (iW). Whitens the
%                  Dirac mode-forward. Bit-identical to BST_INVERSE_LINEAR_2018.
%       STAGE 2 -- Observability SVD (iW*Gm = UL*S*VL') + SNR -> Lambda.
%       STAGE 3 -- Spectral window g(s) = Wiener filter factor on the singular values.
%       STAGE 4 -- Resolution modes W = Phi*VL (one reconstruct); kernel = W*diag(g)*UL'.
%       STAGE 5 -- Source normalization (InverseMeasure) + fold in whitener:
%                  'amplitude' (min-norm) | 'dspm2018' (noise-normalized statistic)
%                  | 'sloreta' (resolution-normalized statistic).
%       TODO -- depth weighting (source covariance WQ); loose orientation;
%               process/panel wrapper.
%
% INPUTS:
%   HeadModel : the standard UNCONSTRAINED surface head model (preferred, same as
%               bst_inverse_linear_2018): .Gain [nCh x 3*nVert], .SurfaceFile.
%               It is transformed to the Dirac eigenbasis internally (bst_dirac),
%               and the vertex leadfield is retained for source-space normalizations.
%               A pre-composed Dirac mode head model (.isDiracEigenmode=1) is also
%               accepted, but then the vertex leadfield is unavailable (sLORETA /
%               depth weighting will not be possible).
%               The Gain rows must already be restricted to the requested channels
%               (same convention as bst_inverse_linear_2018: the caller selects).
%   OPTIONS :
%     .NoiseCovMat   : struct with .NoiseCov [nCh x nCh] (over the requested channels),
%                      and (for 'shrink') .FourthMoment, .nSamples
%     .ChannelTypes  : {1 x nCh} channel type string per Gain row (e.g. 'MEG','MEG GRAD')
%     .NoiseMethod   : {'reg'(default),'median','diag','shrink','none'}
%     .NoiseReg      : ridge fraction for 'reg' (default 0.1, == Brainstorm default)
%     .SnrMethod/.SnrFixed/.SnrRms : SNR -> Lambda regularization (default 'fixed', 3)
%     .InverseMeasure : {'amplitude'(default),'dspm2018','sloreta'}
%     .nModes/.Tau    : Dirac transform params when a vertex head model is given (400, 0.5)
%
% OUTPUTS (STAGE 1):
%   Results :
%     .Whitener       [nCh x nCh]   inverse noise whitener iW_noise (== MNE iW_noise)
%     .GainWhitened   [nCh x nMode] whitened mode-forward iW * Gm
%     .NoiseRankKept  scalar        rank retained in the noise covariance
%     .ChannelTypes, .ModeHemisphere, .Eigenvalues, .HemiGlobalVertices,
%     .DiracEigenFile, .DiracTau, .SurfaceFile  (passed through for later stages)
%
% SEE ALSO: bst_inverse_linear_2018, bst_dirac, tess_eigen
%
% Authors: Diellor Basha, 2026
%
% The two local helpers truncate_and_regularize_covariance and cov1para_local are
% copied verbatim from bst_inverse_linear_2018 (J.C. Mosher et al.) so the whitener
% is bit-identical to the canonical MNE. See that file for full provenance.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@

    % ---- defaults (match bst_inverse_linear_2018) ----
    if (nargin < 2) || isempty(OPTIONS), OPTIONS = struct(); end
    Def = struct('NoiseMethod','reg', 'NoiseReg',0.1, 'ChannelTypes',[], 'NoiseCovMat',[], ...
                 'SnrMethod','fixed', 'SnrFixed',3, 'SnrRms',1000, 'InverseMeasure','amplitude', ...
                 'nModes',400, 'Tau',0.5);
    fn = fieldnames(Def);
    for i = 1:numel(fn)
        if ~isfield(OPTIONS, fn{i}) || isempty(OPTIONS.(fn{i}))
            OPTIONS.(fn{i}) = Def.(fn{i});
        end
    end

    % ---- resolve the head model. The native input is the Dirac mode head model;
    %      a standard surface head model is also accepted and transformed here.
    %      The vertex leadfield is NOT needed: when a source-space normalization
    %      (sLORETA) needs it, the self-consistent EFFECTIVE leadfield is
    %      reconstructed from the modes (bst_dirac Reconstruct of the mode-forward
    %      = the leadfield projected onto the mode span). ----
    if isfield(HeadModel,'isDiracEigenmode') && ~isempty(HeadModel.isDiracEigenmode) && HeadModel.isDiracEigenmode
        CompHM = HeadModel;                       % already a Dirac mode head model
    else
        if ~isfield(HeadModel,'Gain') || isempty(HeadModel.Gain) || mod(size(HeadModel.Gain,2),3) ~= 0
            error('bst_inverse_dirac:notSurface', ...
                'HeadModel must be a Dirac mode head model or an UNCONSTRAINED surface head model.');
        end
        CompHM = bst_dirac(HeadModel, 'nModes', OPTIONS.nModes, 'Tau', OPTIONS.Tau);
    end
    Gm  = double(CompHM.Gain);                    % [nCh x nMode] mode-forward
    nCh = size(Gm, 1);
    if isempty(OPTIONS.NoiseCovMat) || ~isfield(OPTIONS.NoiseCovMat,'NoiseCov') || isempty(OPTIONS.NoiseCovMat.NoiseCov)
        error('bst_inverse_dirac:noNoiseCov', 'OPTIONS.NoiseCovMat.NoiseCov is required.');
    end
    if size(OPTIONS.NoiseCovMat.NoiseCov,1) ~= nCh
        error('bst_inverse_dirac:sizeMismatch', ...
            'NoiseCov is %dx%d but the mode-forward has %d channels.', ...
            size(OPTIONS.NoiseCovMat.NoiseCov,1), size(OPTIONS.NoiseCovMat.NoiseCov,2), nCh);
    end
    if isempty(OPTIONS.ChannelTypes)
        OPTIONS.ChannelTypes = repmat({'MEG'}, 1, nCh);   % single-modality default
    end
    if numel(OPTIONS.ChannelTypes) ~= nCh
        error('bst_inverse_dirac:chanTypes', 'ChannelTypes must have one entry per Gain row (%d).', nCh);
    end

    % =====================================================================
    % STAGE 1 -- NOISE COVARIANCE REGULARIZATION + WHITENER
    %   Faithful port of bst_inverse_linear_2018: zero cross-modality terms,
    %   truncate+regularize each modality separately, assemble iW_noise.
    % =====================================================================
    CROSS_COVARIANCE_CHANNELTYPES = false;   % matches the canonical MNE common path

    C_noise   = OPTIONS.NoiseCovMat.NoiseCov;
    Var_noise = diag(C_noise);

    % Detect a (numerically) diagonal covariance and force 'diag' (MNE behaviour)
    if (norm(C_noise,'fro') - norm(Var_noise,'fro')) < eps(single(norm(Var_noise,'fro')))
        fprintf('BST_INVERSE_DIRAC > Detected diagonal noise covariance, enforcing NoiseMethod="%s".\n', OPTIONS.NoiseMethod);
        OPTIONS.NoiseMethod = 'diag';
    end
    if strcmpi('diag', OPTIONS.NoiseMethod)
        C_noise = diag(diag(C_noise));
    end

    % Channel-type partition
    Unique_ChannelTypes = unique(OPTIONS.ChannelTypes);
    ndx_Channel_Types = cell(1, numel(Unique_ChannelTypes));
    for i = 1:numel(Unique_ChannelTypes)
        ndx_Channel_Types{i} = find(strcmpi(Unique_ChannelTypes(i), OPTIONS.ChannelTypes));
    end

    % Zero cross-modality terms (keep within-modality blocks)
    Cw_noise = zeros(size(C_noise));
    if CROSS_COVARIANCE_CHANNELTYPES
        ndx = [ndx_Channel_Types{:}];
        Cw_noise(ndx,ndx) = C_noise(ndx,ndx);
    else
        for i = 1:numel(Unique_ChannelTypes)
            ndx = ndx_Channel_Types{i};
            Cw_noise(ndx,ndx) = C_noise(ndx,ndx);
        end
    end
    Cw_noise = (Cw_noise + Cw_noise')/2;

    % FourthMoment / nSamples for 'shrink'
    FourthMoment = zeros(size(C_noise));
    nSamples = [];
    if strcmpi(OPTIONS.NoiseMethod, 'shrink')
        if ~isfield(OPTIONS.NoiseCovMat,'FourthMoment') || isempty(OPTIONS.NoiseCovMat.FourthMoment)
            error('bst_inverse_dirac:shrink', 'NoiseMethod ''shrink'' requires NoiseCovMat.FourthMoment.');
        end
        FourthMoment = OPTIONS.NoiseCovMat.FourthMoment;
        if ~isfield(OPTIONS.NoiseCovMat,'nSamples') || isempty(OPTIONS.NoiseCovMat.nSamples)
            error('bst_inverse_dirac:shrink', 'NoiseMethod ''shrink'' requires NoiseCovMat.nSamples.');
        end
        nSamples = OPTIONS.NoiseCovMat.nSamples(1);
    end

    % Truncate + regularize each modality, assemble the inverse whitener
    iWw_noise = zeros(size(Cw_noise));
    NoiseRankKept = 0;
    for i = 1:numel(Unique_ChannelTypes)
        ndx = ndx_Channel_Types{i};
        [Cw_noise(ndx,ndx), iWw_noise(ndx,ndx)] = truncate_and_regularize_covariance( ...
            Cw_noise(ndx,ndx), OPTIONS.NoiseMethod, Unique_ChannelTypes{i}, ...
            OPTIONS.NoiseReg, FourthMoment(ndx,ndx), nSamples);
        NoiseRankKept = NoiseRankKept + sum(abs(diag(iWw_noise(ndx,ndx))) > 0);
    end
    iW_noise = iWw_noise;

    % ---- whiten the Dirac mode-forward ----
    GainWhitened = iW_noise * Gm;

    fprintf('BST_INVERSE_DIRAC > STAGE 1 done: noise rank kept %d/%d, whitened mode-forward [%d x %d].\n', ...
        NoiseRankKept, nCh, size(GainWhitened,1), size(GainWhitened,2));

    % =====================================================================
    % STAGE 2 -- OBSERVABILITY SVD + SNR -> LAMBDA
    %   One SVD of the whitened mode-forward:  Lt = iW*Gm = UL * S * VL'.
    %   The singular values SL are the OBSERVABILITY spectrum (how strongly the
    %   sensors see each independent source direction); UL (sensor) and VL (mode)
    %   are the resolution axes. Stages 3-5 are all DIAGONAL windows on this one
    %   spectrum -- the filter-factor view of linear inverse methods.
    %   Source prior is the unit prior (WQ = I); Lambda sets the SNR corner.
    % =====================================================================
    [UL, Ssvd, VL] = svd(GainWhitened, 'econ');   % [nCh x r], [r x r], [nMode x r]
    SL  = diag(Ssvd);                             % observability singular values
    SL2 = SL.^2;
    Rank_Leadfield = sum(SL > length(SL)*eps(single(SL(1))));

    switch lower(OPTIONS.SnrMethod)
        case 'rms'
            Lambda = OPTIONS.SnrRms^2;  SNR = Lambda * SL2(1);
        case 'fixed'
            SNR = OPTIONS.SnrFixed^2;   Lambda = SNR / mean(SL2);   % Hamalainen mean-eigenvalue
        otherwise
            error('bst_inverse_dirac:snr', 'Unsupported SnrMethod="%s".', OPTIONS.SnrMethod);
    end

    fprintf('BST_INVERSE_DIRAC > STAGE 2 done: whitened-leadfield rank %d/%d, assumed SNR=%.1f (%.1f dB), Lambda=%.3e.\n', ...
        Rank_Leadfield, length(SL), SNR, 10*log10(SNR), Lambda);

    % =====================================================================
    % STAGE 3 -- SPECTRAL WINDOW (filter factor) on the observability axis
    %   Every linear inverse is a window g(s) applied to each singular direction.
    %   Amplitude min-norm is the Wiener/Tikhonov window:
    %       g(s) = Lambda*s / (Lambda*s^2 + 1)
    %   The mode-coefficient kernel is then VL*diag(g)*UL' (kept factored).
    % =====================================================================
    gWin = (Lambda * SL) ./ (Lambda * SL2 + 1);   % [r x 1] Wiener window on s
    KernelMode = VL * (gWin .* UL');              % [nMode x nCh] mode kernel (whitened data)

    fprintf('BST_INVERSE_DIRAC > STAGE 3 done: Wiener window g(s), half-power at i=%d/%d.\n', ...
        sum(Lambda*SL2 > 1), numel(SL));

    % =====================================================================
    % STAGE 4 -- VERTEX-SPACE RESOLUTION MODES  W = Phi * VL
    %   Reconstruct the r right singular vectors ONCE (bst_dirac Reconstruct,
    %   J = Phi*c). The whitened-data current kernel is then the rank-r product
    %       KvtxW = W * diag(g) * UL'                 [3*nVert x nCh]
    %   ( = Phi-expansion of the mode kernel VL*diag(g)*UL' ). Because Phi is
    %   B-orthonormal the underlying min-norm minimizes ||c||_2 = ||J||_B (the
    %   mass/area-weighted continuous min-norm). iW is folded in at STAGE 5.
    % =====================================================================
    Wres  = bst_dirac(CompHM, 'Reconstruct', VL')';   % [3*nVert x r] vertex resolution modes
    nVert = size(Wres,1) / 3;
    KvtxW = Wres * (gWin .* UL');                      % [3*nVert x nCh] whitened-data kernel

    fprintf('BST_INVERSE_DIRAC > STAGE 4 done: %d resolution modes, current kernel [%d x %d] (whitened-data).\n', ...
        size(Wres,2), size(KvtxW,1), size(KvtxW,2));

    % =====================================================================
    % STAGE 5 -- SOURCE NORMALIZATION (measure) + fold in the whitener
    %   Per-source (unconstrained 3-component) normalization, computed on the
    %   WHITENED-data vertex kernel (so the source-noise covariance is the kernel
    %   row gram). Port of bst_inverse_linear_2018 dspm2018 (lines 862-901).
    %     amplitude : no normalization (min-norm current)
    %     dspm2018  : divide each source by its projected noise std
    %                 sqrt(trace(K_v K_v')) = noise-normalized statistic (depth-bias free)
    %     sloreta   : normalize by the resolution R_v^{-1/2} (PM 2002; needs Gvtx)
    % =====================================================================
    switch lower(OPTIONS.InverseMeasure)
        case 'amplitude'
            Knorm = KvtxW;                                     % min-norm current (no normalization)
        case 'dspm2018'
            rowN2 = sum(KvtxW.^2, 2);                          % [3*nVert x 1] row norm^2 (white noise)
            Snorm = sqrt(sum(reshape(rowN2, 3, []), 1))';      % [nVert x 1] sqrt(trace per source)
            Snorm3 = reshape(repmat(Snorm', 3, 1), [], 1); Snorm3(Snorm3 < eps) = eps;
            Knorm = KvtxW ./ Snorm3;
        case 'sloreta'
            % per-source resolution R_v = K_v * L_v (3x3); normalize by R_v^{-1/2}.
            % The effective leadfield is RECONSTRUCTED from the modes (the forward
            % the inverse actually inverts), so no vertex head model is needed.
            Lw = iW_noise * bst_dirac(CompHM, 'Reconstruct', Gm);   % [nCh x 3*nVert] whitened
            Knorm = zeros(size(KvtxW));
            for v = 1:nVert
                idx = (v-1)*3 + (1:3);
                Rv  = KvtxW(idx,:) * Lw(:,idx);                % 3x3 resolution
                [Ur,Sr,Vr] = svd((Rv+Rv')/2);  sr = diag(Sr);
                rk = sum(sr > length(sr)*eps(single(sr(1))));
                SIR = Vr(:,1:rk) * diag(1./sqrt(sr(1:rk))) * Ur(:,1:rk)';   % R_v^{-1/2}
                Knorm(idx,:) = SIR * KvtxW(idx,:);
            end
        otherwise
            error('bst_inverse_dirac:measure', ...
                'InverseMeasure must be ''amplitude'', ''dspm2018'' or ''sloreta'', got ''%s''.', ...
                OPTIONS.InverseMeasure);
    end
    ImagingKernel = Knorm * iW_noise;                         % [3*nVert x nCh] on RAW data

    fprintf('BST_INVERSE_DIRAC > STAGE 5 done: %s imaging kernel [%d x %d].\n', ...
        OPTIONS.InverseMeasure, size(ImagingKernel,1), size(ImagingKernel,2));

    % =====================================================================
    % Assemble result (carry mode metadata + stage products for later stages)
    % =====================================================================
    Results = struct();
    Results.Whitener      = iW_noise;          % STAGE 1
    Results.GainWhitened  = GainWhitened;      % STAGE 1
    Results.NoiseRankKept = NoiseRankKept;     % STAGE 1
    Results.NoiseMethod   = OPTIONS.NoiseMethod;
    Results.ChannelTypes  = OPTIONS.ChannelTypes;
    Results.UL            = UL;                % STAGE 2: left singular vectors (sensor space)
    Results.SL            = SL;                % STAGE 2: whitened-leadfield singular values
    Results.Lambda        = Lambda;           % STAGE 2: source-variance prior (regularizer)
    Results.SNR           = SNR;              % STAGE 2: assumed SNR (squared)
    Results.RankLeadfield = Rank_Leadfield;   % STAGE 2
    Results.ImagingKernelMode = KernelMode * iW_noise;  % STAGE 3: c_hat = . * d_raw  [nMode x nCh] (amplitude)
    Results.ImagingKernel = ImagingKernel;    % STAGE 5: J = ImagingKernel * d_raw  [3*nVert x nCh]
    Results.InverseMeasure = OPTIONS.InverseMeasure;
    % pass-through Dirac metadata
    for f = {'Eigenvalues','ModeHemisphere','HemiGlobalVertices','DiracEigenFile','DiracTau','SurfaceFile','nModes'}
        if isfield(CompHM, f{1}), Results.(f{1}) = CompHM.(f{1}); end
    end
end


% =============================================================================
% ===== LOCAL HELPERS (verbatim from bst_inverse_linear_2018, J.C. Mosher) =====
% =============================================================================
function [Cov,iW] = truncate_and_regularize_covariance(Cov,Method,Type,NoiseReg,FourthMoment,nSamples)
% Cov is the covariance matrix, to be regularized using Method
% Type is the sensor type for display purposes
% NoiseReg is the regularization fraction, if Method "reg" selected
% FourthMoment and nSamples are used if Method "shrinkage" selected

VERBOSE = true;
Cov = (Cov + Cov')/2;
[Un,Sn2] = svd(Cov,'econ');
Sn = sqrt(diag(Sn2));
tol = length(Sn) * eps(single(Sn(1)));
Rank_Noise = sum(Sn > tol);
if VERBOSE
    fprintf('BST_INVERSE_DIRAC > Rank of the ''%s'' channels, keeping %.0f noise eigenvalues out of %.0f original set\n',...
        Type,Rank_Noise,length(Sn));
end
Un = Un(:,1:Rank_Noise);
Sn = Sn(1:Rank_Noise);
Cov = Un*diag(Sn.^2)*Un';

if VERBOSE
    fprintf('BST_INVERSE_DIRAC > Using the ''%s'' method of covariance regularization.\n',Method);
end

switch(Method) % {'shrink', 'reg', 'diag', 'none', 'median'}
    case 'none'
        iW = Un*diag(1./Sn)*Un';
        if VERBOSE, fprintf('BST_INVERSE_DIRAC > No regularization applied to covariance matrix.\n'); end
    case 'median'
        if VERBOSE
            fprintf('BST_INVERSE_DIRAC > Covariance regularized by flattening tail of eigenvalues spectrum to the median value of %.1e\n',median(Sn));
        end
        Sn = max(Sn,median(Sn));
        Cov = Un*diag(Sn.^2)*Un';
        iW = Un*diag(1./Sn)*Un';
    case 'diag'
        Cov = diag(diag(Cov));
        iW = diag(1./sqrt(diag(Cov)));
        if VERBOSE, fprintf('BST_INVERSE_DIRAC > Covariance matrix reduced to diagonal.\n'); end
    case 'reg'
        RidgeFactor = mean(diag(Sn2)) * NoiseReg;   % Hamalainen mean-eigenvalue ridge
        Cov = Cov + RidgeFactor * eye(size(Cov,1));
        iW = Un*diag(1./sqrt(Sn.^2 + RidgeFactor))*Un';
        if VERBOSE
            fprintf('BST_INVERSE_DIRAC > Diagonal of %.1f%% of the average eigenvalue added to covariance matrix.\n',NoiseReg * 100);
        end
    case 'shrink'
        [Cov,shrinkage]=cov1para_local(Cov,FourthMoment,nSamples);
        if VERBOSE, fprintf('\nShrinkage factor is %f\n\n',shrinkage); end
        [Un,Sn2] = svd(Cov,'econ');
        Sn = sqrt(diag(Sn2));
        tol = length(Sn) * eps(single(Sn(1)));
        Rank_Noise = sum(Sn > tol);
        if VERBOSE
            fprintf('BST_INVERSE_DIRAC > Ledoit covariance regularization, after shrinkage, rank of the %s channels, keeping %.0f noise eigenvalues out of %.0f original set\n',...
                Type,Rank_Noise,length(Sn));
        end
        Un = Un(:,1:Rank_Noise);
        Sn = Sn(1:Rank_Noise);
        Cov = Un*diag(Sn.^2)*Un';
        iW = Un*diag(1./Sn)*Un';
    otherwise
        error(['BST_INVERSE_DIRAC > Unknown covariance regularization method: NoiseMethod="' Method '"']);
end
end

% ======== Ledoit Shrinkage (verbatim from bst_inverse_linear_2018) ========
function [sNoiseCov,shrinkage]=cov1para_local(NoiseCov,FourthMoment,nSamples)
n=size(NoiseCov,1);
meanvar=mean(diag(NoiseCov));
prior=meanvar*eye(n);
phiMat = FourthMoment - NoiseCov.^2;
phi=sum(sum(phiMat));
gamma=norm(NoiseCov-prior,'fro')^2;
kappa=phi/gamma;
shrinkage=max(0,min(1,kappa/nSamples));
sNoiseCov=shrinkage*prior+(1-shrinkage)*NoiseCov;
end
