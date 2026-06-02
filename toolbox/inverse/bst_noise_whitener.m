function iW = bst_noise_whitener(NoiseCov, ChannelTypes, NoiseMethod, NoiseReg, FourthMoment, nSamples)
% BST_NOISE_WHITENER: Per-modality noise whitener iW = C^(-1/2).
%
% USAGE:  iW = bst_noise_whitener(NoiseCov, ChannelTypes, NoiseMethod, NoiseReg, FourthMoment, nSamples)
%
% DESCRIPTION:
%   Verbatim extraction of the per-modality whitener used by
%   bst_inverse_linear_2018 (lines 237-446 + subfunctions
%   truncate_and_regularize_covariance and cov1para_local). Regularizes and
%   whitens each modality (channel type) separately, with cross-modality
%   covariances zeroed, exactly as the standard 2018 inverse does. Kept as a
%   standalone function so the eigenmode inverse can reuse the same whitener
%   without depending on the standard solver.
%
%   NoiseCov     : [nCh x nCh] noise covariance, good channels only
%   ChannelTypes : 1 x nCh cell of channel type strings (e.g. {'MEG MAG',...})
%   NoiseMethod  : 'reg' | 'diag' | 'none' | 'shrink' | 'median'
%   NoiseReg     : scalar in [0,1] (used by 'reg')
%   FourthMoment : [nCh x nCh] (only required for 'shrink'); default zeros
%   nSamples     : scalar (only required for 'shrink'); default []
%
% NOTE: This duplicates logic in bst_inverse_linear_2018. Do not edit the math
%       here independently; if the standard whitener changes, re-sync verbatim.
%       Future consolidation (make bst_inverse_linear_2018 call this) is tracked
%       as tech debt in dev/2026-06-02-eigenmode-interactive-inverse-design.md.
%
% Authors: Diellor Basha, 2026

if (nargin < 5) || isempty(FourthMoment)
    FourthMoment = zeros(size(NoiseCov));
end
if (nargin < 6)
    nSamples = [];
end
% Assemble the OPTIONS-shaped struct the pasted block expects
OPTIONS.NoiseCovMat.NoiseCov     = NoiseCov;
OPTIONS.NoiseCovMat.FourthMoment = FourthMoment;
OPTIONS.NoiseCovMat.nSamples     = nSamples;
OPTIONS.ChannelTypes             = ChannelTypes;
OPTIONS.NoiseMethod              = NoiseMethod;
OPTIONS.NoiseReg                 = NoiseReg;

% Verbatim from bst_inverse_linear_2018.m line 196
CROSS_COVARIANCE_CHANNELTYPES = false;

% ===== BEGIN VERBATIM COPY: bst_inverse_linear_2018.m lines 237-446 =====
C_noise = OPTIONS.NoiseCovMat.NoiseCov; % all of the channels requested
Var_noise = diag(C_noise); % Diagonal vector of the noise variances per channel
nChannels = length(Var_noise);

% JCM: March 2015, probably don't need this, but will retain
% Detect if the input noise covariance matrix is or should be diagonal
if (norm(C_noise,'fro') - norm(Var_noise,'fro')) < eps(single(norm(Var_noise,'fro')))
    % no difference between the full matrix and the diagonal matrix
    disp(['BST_INVERSE > Detected diagonal noise covariance, enforcing option NoiseMethod="' OPTIONS.NoiseMethod '".']);
    OPTIONS.NoiseMethod = 'diag';
end

if strcmpi('diag',OPTIONS.NoiseMethod)
    C_noise = diag(diag(C_noise)); % force matrix to be diagonal, that's all the user wants
end

% Commentary April 2015: There is a tendency to apply generic scale values to the
% channels to balance them with regards to units, such as Volts and Tesla.
% But we also have a problem with gradiometer channels vs magnetometers. So
% the "natural" way is to use the channel variances themselves to balance
% out the differences between modalities. But we don't want to do each
% channel, since a dead channel has (near) zero variance. So instead we
% calculate a common variance for each modality, to get us in the ball
% park. So we initially treat each modality as Independent and Identically
% Distributed (IID) and pre-whiten by this to bring the modalities into
% closer alignment.

% Because the units can be different, we need to first balance the
% different types of arrays. How many unique array types do we have?

% Updated Commentary August 2016, by John Mosher: The issue of balancing is
% primarily a problem of multi-modal statistics, i.e., across multiple
% channel types. Ideally, we should allow for possible cross-covariances
% between the modalities to assist in the overal spatial correlation of the
% array. The problem is that we simultaneously need to "regularize" the
% covariance matrices to catch bad channels and experimental dependencies
% that creep into measurements, while also calculating cross-dependencies
% between modalities of disparate units.
%
% We have tried numerous "black-box" methods to get, for example, Neuromag
% magnetometers and gradiometers combined in the same covariance matrix.
% Unfortunately, each situation can be unique, such as in a clinical
% setting, where the magnetometers can span an enormous dynamic range,
% while gradiometers span a much tighter range, and these ranges overlap.
%
% In consultation with co-investigator Matti Hamalainen and MNE, we opted
% for now to keep the use of the multimodal calculation simpler by
% eliminating the cross modality terms, then regularizing within each
% modality. This has been the default approach to multimodal noise
% covariance calculations in Brainstorm since 2011, and we retain that
% here.

% So how many channel types do we have:
Unique_ChannelTypes = unique(OPTIONS.ChannelTypes);

% What are their indices
ndx_Channel_Types = cell(1,length(Unique_ChannelTypes)); % index for each channel type
for i = 1:length(Unique_ChannelTypes)
    ndx_Channel_Types{i} = find(strcmpi(Unique_ChannelTypes(i),OPTIONS.ChannelTypes));
end

% not sure if Ledoit's shrinkage method will work across mixed modalities
if length(Unique_ChannelTypes) > 1 && strcmpi('shrink',OPTIONS.NoiseMethod) && CROSS_COVARIANCE_CHANNELTYPES,
    fprintf('NOTE: Noise Regularization ''shrink'' selected with multiple channel types, not sure that will work with cross covariances.\n');
end


% Initialize the modality whitening matrices:
% iPrior_Matrix = eye(nChannels);

% Retain the below code in comments if only to discuss what does and does
% not work in trying to capture cross-modalities. August 2016

%     % calculate the average variance per channel as a prior
%     Prior_IID = zeros(length(Unique_ChannelTypes),1);
%     Prior_Vector = zeros(nChannels,1); % initialize as a vector
%
%     for i = 1:length(Unique_ChannelTypes),
%         ndx = find(strcmpi(Unique_ChannelTypes(i),OPTIONS.ChannelTypes)); % which ones
%         Prior_IID(i) = mean(Var_noise(ndx)); % mean of this type
%         % let's get a bit more sophisticated on calculating this IID value
%         % We wouldn't want a few extreme values distorting our mean
%         % How many channels are there:
%         len_ndx = length(ndx);
%         % let's toss out the upper and lower values
%         ndx_clip = round(len_ndx/10); % 10%
%         Variances_This_Modality = sort(Var_noise(ndx));
%         % trim to central values
%         Variances_This_Modality = Variances_This_Modality((ndx_clip+1):(end-ndx_clip));
%         Prior_IID(i) = median(Variances_This_Modality); % mean of the middle values
%         Prior_Vector(ndx) = sqrt(Prior_IID(i)); % map to this part of the array
%     end
%
%     %TODO: Test for bizarre case of Prior_IID too small
%
%     % build whitener to balance out the different types of channels
%     iPrior_Matrix = spdiags(1./Prior_Vector,0,nChannels,nChannels);
%
%     % Now we can use this iPrior_Matrix to "pre-whiten" the noise covariance
%
%     % Block Whitened noise covariance
%     Cw_noise = iPrior_Matrix * C_noise * iPrior_Matrix;
%
%     % Now the units imbalance between different subarrays is theoretically
%     % balanced.

% January 2018, for consistency with Hamalainen's MNE and with older
% Brainstorm code, Mosher making some changes to how noise covariance is
% calculated and regularized.

% Now remove the off_diagonal components between modalities, if desired
Cw_noise = zeros(size(C_noise)); % initialize

if CROSS_COVARIANCE_CHANNELTYPES % do allow for cross covariances

    ndx = [ndx_Channel_Types{:}]; % all of the channels
    Cw_noise(ndx,ndx) = C_noise(ndx,ndx); % including cross terms

else % don't allow for cross covariances

    for i = 1:length(Unique_ChannelTypes)
        % for each unique modality
        ndx = ndx_Channel_Types{i}; % which ones
        Cw_noise(ndx,ndx) = C_noise(ndx,ndx);
    end
    % the off diagonal components corresponding to cross-modality
    % covariances have been removed

end

% Mosher, Feb 2018 change, we may want the cross modalities, in order to
% exploit any additional information between them. What we should have is a
% "user flag", but for now I will hard-wire it in and discuss with the
% others. Matti uses the cross information between GRADS and MAGS, but not
% between the other modalities. However, that is the point of the advance
% e-physiological studies, merging modalities.

% By continuing to call it "Cw_noise", the below codes can continue as
% before.

% Before any additional regularizations that may interfere

% Make sure the noise covariance matrix is strictly symmetric
% (the previous operations may cause rounding errors that make the matrix not exactly symmetrical)
Cw_noise = (Cw_noise + Cw_noise')/2;

% So the block whitened noise covariance matrix is Cw_noise.

% If the modalities were truly balanced, we could proceed with
% regularization in one fell swoop; however, the reality is that
% regularization should be applied modality by modality, with the cross
% covariance of the modalities zeroed out.

% initialize inverse whitener for noise covariance matrix
iWw_noise = zeros(size(Cw_noise));

FourthMoment = zeros(size(C_noise)); % initialize
nSamples = [];
if strcmpi(OPTIONS.NoiseMethod, 'shrink')
   % Has the user calculated the Fourth Order moments?
   if ~isfield(OPTIONS.NoiseCovMat,'FourthMoment')
      error('BST_INVERSE > For Method ''shrink'' please recalculate Noise Covariance, to include Fourth Order Moments');
   else
      FourthMoment = OPTIONS.NoiseCovMat.FourthMoment;
   end

   % How many samples used to calculate covariance
   if isempty(OPTIONS.NoiseCovMat.nSamples)
      error('BST_INVERSE > No noise samples found. For Method ''shrink'' please recalculate Noise Covariance from actual data samples.');
   else
      nSamples = OPTIONS.NoiseCovMat.nSamples(1);
      % FIX: What if different lengths are used, but only an issue at this
      % point for "shrink" method
   end
end

% Mosher, Feb 2018, now the trick to regularization across modalities is to
% apply the regularization within the modality itself, but not across
% modalites. In other words, calculate the average variance of the MAGS
% separate from the GRADS separate from SEEG.

% First, truncate and regularize each modality separately
for i = 1:length(Unique_ChannelTypes)
   % each modality
    ndx = ndx_Channel_Types{i}; % which ones

    % regularize and form whitener for each modality, using local
    % subfunction
    [Cw_noise(ndx,ndx),iWw_noise(ndx,ndx)] = ...
        truncate_and_regularize_covariance(Cw_noise(ndx,ndx),...
        OPTIONS.NoiseMethod,Unique_ChannelTypes{i},...
        OPTIONS.NoiseReg,FourthMoment(ndx,ndx),nSamples);
end

% So now we have calculated the noise covariance matrix across all
% modalities, but have eliminated the cross covariances between modalities.
% Within each modality, we have checked for rank deficiencies and have
% added possible regularization within the modality.

% now, if cross modality is desired, we need to calculate the overall
% inverse, using the regularized submatrices (all which may have had their
% diagonal terms altered by regularization), but now we don't apply any
% additional regularization to the overall covariance matrix

if length(Unique_ChannelTypes) > 1 && CROSS_COVARIANCE_CHANNELTYPES % we do want the cross terms in the inverse
    ndx = [ndx_Channel_Types{:}]; % all of the channels
    % don't regularize, may still need truncation if deficient
    [Cw_noise(ndx,ndx),iWw_noise(ndx,ndx)] = ...
        truncate_and_regularize_covariance(Cw_noise(ndx,ndx),'none','ALL');
end
% ===== END VERBATIM COPY =====

iW = iWw_noise;
end

% ===== BEGIN VERBATIM COPY: bst_inverse_linear_2018.m lines 1098-1232 =====
function [Cov,iW] = truncate_and_regularize_covariance(Cov,Method,Type,NoiseReg,FourthMoment,nSamples)
% Cov is the covariance matrix, to be regularized using Method
% Type is the sensor type for display purposes
% NoiseReg is the regularization fraction, if Method "reg" selected
% FourthMoment and nSamples are used if Method "shrinkage" selected

VERBOSE = true; % be talkative about what's happening

% Ensure symmetry
Cov = (Cov + Cov')/2;

% Note,impossible to be complex by above symmetry check
% Decompose just this covariance.
[Un,Sn2] = svd(Cov,'econ');
Sn = sqrt(diag(Sn2)); % singular values
tol = length(Sn) * eps(single(Sn(1))); % single precision tolerance
Rank_Noise = sum(Sn > tol);

if VERBOSE
    fprintf('BST_INVERSE > Rank of the ''%s'' channels, keeping %.0f noise eigenvalues out of %.0f original set\n',...
        Type,Rank_Noise,length(Sn));
end

Un = Un(:,1:Rank_Noise);
Sn = Sn(1:Rank_Noise);

% now rebuild the noise covariance matrix with just the non-zero
% components
Cov = Un*diag(Sn.^2)*Un'; % possibly deficient matrix now

% With this modality truncated, see if we need any additional
% regularizations, and build the inverse whitener

if VERBOSE
    fprintf('BST_INVERSE > Using the ''%s'' method of covariance regularization.\n',Method);
end

switch(Method) % {'shrink', 'reg', 'diag', 'none', 'median'}

    case 'none'
        %  "none" in Regularization means no
        % regularization was applied to the computed Noise Covariance
        % Matrix
        % Do Nothing to Cw_noise
        iW = Un*diag(1./Sn)*Un'; % inverse whitener
        if VERBOSE
            fprintf('BST_INVERSE > No regularization applied to covariance matrix.\n');
        end


    case 'median'
        if VERBOSE
            fprintf('BST_INVERSE > Covariance regularized by flattening tail of eigenvalues spectrum to the median value of %.1e\n',median(Sn));
        end
        Sn = max(Sn,median(Sn)); % removes deficient small values
        Cov = Un*diag(Sn.^2)*Un'; % rebuild again.
        iW = Un*diag(1./Sn)*Un'; % inverse whitener

    case 'diag'
        Cov = diag(diag(Cov)); % strip to diagonal
        iW = diag(1./sqrt(diag(Cov))); % inverse of diagonal
        if VERBOSE
            fprintf('BST_INVERSE > Covariance matrix reduced to diagonal.\n');
        end

    case 'reg'
        % The unit of "Regularize Noise Covariance" is as a percentage of
        % the mean variance of the modality.

        % Ridge Regression:
        % Commented out Feb 2018 in favor of the mean eigenvalue, to align
        % with Hamalainen.
        % RidgeFactor = Sn2(1) * NoiseReg ; % percentage of max.
        % Use instead this one:
        RidgeFactor = mean(diag(Sn2)) * NoiseReg; % Hamalainen's preferred measure
        %(Note, the mean of the eigenvalues is the mean of the diagonal
        %values).

        Cov = Cov + RidgeFactor * eye(size(Cov,1));
        % wrong: iW = Un*diag(1./(Sn + sqrt(RidgeFactor)))*Un'; % inverse whitener
        % Fixed Feb 2018:
        iW = Un*diag(1./sqrt(Sn.^2 + RidgeFactor))*Un'; % inverse whitener, symmetric

        if VERBOSE
            fprintf('BST_INVERSE > Diagonal of %.1f%% of the average eigenvalue added to covariance matrix.\n',NoiseReg * 100);
        end


    case 'shrink'
        % Method of Ledoit, recommended by Alexandre Gramfort

        % Need to scale the Fourth Moment for the modalities

        % use modified version of cov1para attached to this function
        % TODO, can we adjust this routine to handle different numbers of
        % samples in the generation of the fourth order moments
        % calculation? As of August 2016, still relying on a single scalar
        % number.
        [Cov,shrinkage]=cov1para_local(Cov,FourthMoment,nSamples);
        if VERBOSE
            fprintf('\nShrinkage factor is %f\n\n',shrinkage)
        end
        % we now have the "shrunk" whitened noise covariance
        % Recalculate
        [Un,Sn2] = svd(Cov,'econ');
        Sn = sqrt(diag(Sn2)); % singular values
        tol = length(Sn) * eps(single(Sn(1))); % single precision tolerance
        Rank_Noise = sum(Sn > tol);

        if VERBOSE
            fprintf('BST_INVERSE > Ledoit covariance regularization, after shrinkage, rank of the %s channels, keeping %.0f noise eigenvalues out of %.0f original set\n',...
                Type,Rank_Noise,length(Sn));
        end

        Un = Un(:,1:Rank_Noise);
        Sn = Sn(1:Rank_Noise);

        % now rebuild the noise covariance matrix with just the non-zero
        % components
        Cov = Un*diag(Sn.^2)*Un'; % possibly deficient matrix now

        iW = Un*diag(1./Sn)*Un'; % inverse whitener


    otherwise
        error(['BST_INVERSE > Unknown covariance regularization method: NoiseMethod="' Method '"']);

end % method of regularization


% Note the design of full rotating whiteners. We don't expect dramatic reductions in
% rank here, and it's convenient to rotate back to the original space.
% Note that these whitener matrices may not be of full rank.

end
% ===== END VERBATIM COPY =====

% ===== BEGIN VERBATIM COPY: bst_inverse_linear_2018.m lines 1236-1361 =====
function [sNoiseCov,shrinkage]=cov1para_local(NoiseCov,FourthMoment,nSamples)
% Based on Ledoit's "cov1para" with some modifications
%   x is t x n, returns
%   sigma n x n
%
% shrinkage is the final computed shrinkage factor, used to weight the
%  i.i.d. prior vs the sample estimate. If shrinkage is specified, it is
%  used on input; else, it's computed.

% Original code from
% http://www.ledoit.net/cov1para.m
% Original Ledoit comments:
% function sigma=cov1para(x)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards one-parameter matrix:
%    all variances are the same
%    all covariances are zero
% if shrink is specified, then this const. is used for shrinkage

% Based on
% http://www.ledoit.net/ole1_abstract.htm
% http://www.ledoit.net/ole1a.pdf (PDF of paper)
%
% A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices
% Olivier Ledoit and Michael Wolf
% Journal of Multivariate Analysis, Volume 88, Issue 2, February 2004, pages 365-411
%
% Abstract
% Many economic problems require a covariance matrix estimator that is not
% only invertible, but also well-conditioned (that is, inverting it does
% not amplify estimation error). For large-dimensional covariance matrices,
% the usual estimator - the sample covariance matrix - is typically not
% well-conditioned and may not even be invertible. This paper introduces an
% estimator that is both well-conditioned and more accurate than the sample
% covariance matrix asymptotically. This estimator is distribution-free and
% has a simple explicit formula that is easy to compute and interpret. It
% is the asymptotically optimal convex combination of the sample covariance
% matrix with the identity matrix. Optimality is meant with respect to a
% quadratic loss function, asymptotically as the number of observations and
% the number of variables go to infinity together. Extensive Monte-Carlo
% confirm that the asymptotic results tend to hold well in finite sample.


% Original Code Header, updated to be now from (2014)
% http://www.econ.uzh.ch/faculty/wolf/publications/cov1para.m.zip
%
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards constant correlation matrix
% if shrink is specified, then this constant is used for shrinkage
%
% The notation follows Ledoit and Wolf (2003, 2004)
% This version 04/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is released under the BSD 2-clause license.
%
% Copyright (c) 2014, Olivier Ledoit and Michael Wolf
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wolf's site now,
% http://www.econ.uzh.ch/faculty/wolf/publications/cov1para.m.zip
% some differences from original Ledoit that confirm Mosher's
% original re-coding.

% % de-mean returns
% [t,n]=size(x);
% meanx=mean(x);
% x=x-meanx(ones(t,1),:);

% compute sample covariance matrix
% Provided
% NoiseCov=(1/t).*(x'*x);

% compute prior
n=size(NoiseCov,1); % number of channels
meanvar=mean(diag(NoiseCov));
prior=meanvar*eye(n); % Note, should be near identity by our pre-whitening

% what we call p
%y=x.^2;
%phiMat=y'*y/t - NoiseCov.^2;
phiMat = FourthMoment - NoiseCov.^2;
phi=sum(sum(phiMat));

% what we call r is not needed for this shrinkage target

% what we call c
gamma=norm(NoiseCov-prior,'fro')^2;

% compute shrinkage constant
kappa=phi/gamma;
% ensure bounded between zero and one
shrinkage=max(0,min(1,kappa/nSamples));

% compute shrinkage estimator
sNoiseCov=shrinkage*prior+(1-shrinkage)*NoiseCov;

end
% ===== END VERBATIM COPY =====
