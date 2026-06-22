function [Kernel, Info] = bst_eigenmodes_transform(Gain, Phi, varargin)
% BST_EIGENMODES_TRANSFORM: Unregularized sensor->eigenmode spatial transform.
%
% USAGE:  [Kernel, Info] = bst_eigenmodes_transform(Gain, Phi)
%         [Kernel, Info] = bst_eigenmodes_transform(Gain, Phi, 'Tol', tol)
%
% DESCRIPTION:
%     Builds the composite transform that maps sensor recordings directly to
%     Laplace-Beltrami eigenmode coefficients, with NO regularization. This is
%     the spatial analogue of a Fourier transform: a fixed change of basis with
%     no tuning parameters. Regularization/denoising is a separate, optional
%     step applied to the coefficients afterwards (a future filter library).
%
%     The compressed lead field is L_tilde = Gain * Phi, where each column is
%     the sensor topography of one eigenmode. The transform is the Moore-Penrose
%     pseudoinverse computed via SVD (L_tilde = U*S*V'):
%
%         Kernel = pinv(L_tilde) = V * diag(1./s) * U'      [K x nch]
%         Theta  = Kernel * Data                            [K x nTime]
%
%     SVD is used (rather than the normal equations) because the correct closed
%     form depends on the K/nch ratio -- the left-inverse (K<=nch) and the
%     right-inverse (K>=nch) each have a singular normal matrix in the other
%     regime -- and because forming Gain*Phi's Gram matrix would square an
%     already-large condition number (high-lambda modes are nearly invisible to
%     the sensors). Small singular values are floored (rank guard) but not
%     otherwise weighted: the transform is unregularized but rank-safe.
%
% INPUTS:
%     Gain : [nch x nVert] constrained (fixed-orientation) lead field, already
%            restricted to the channels to use.
%     Phi  : [nVert x K]   eigenmode matrix (caller truncates to K modes).
%
% OPTIONS (name-value):
%     'Tol' : singular-value floor (rank guard). Default [] -> MATLAB pinv
%             default max(size(L_tilde))*eps(max(s)).
%
% OUTPUTS:
%     Kernel : [K x nch] transform A = pinv(L_tilde).
%     Info   : struct with fields
%              .CompressedLF    [nch x K]  L_tilde = Gain*Phi
%              .SingularValues  [min(nch,K) x 1]  all singular values of L_tilde (descending)
%              .Rank            scalar     number of singular values > Tol
%              .ConditionNumber scalar     s(1)/s(Rank)
%              .Tol             scalar     floor used
%              .nModes          scalar     K
%
% SEE ALSO: bst_eigenmodes_project, bst_inverse_eigenmodes, in_tess_eigenmodes

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
%
% Authors: Diellor Basha, 2026

%% ===== PARSE OPTIONS =====
Tol = [];
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'tol'
            Tol = varargin{i+1};
            if ~isempty(Tol) && (~isscalar(Tol) || Tol < 0)
                error('bst_eigenmodes_transform: Tol must be a non-negative scalar.');
            end
    end
end

Gain = double(Gain);
Phi  = double(Phi);

%% ===== COMPRESS LEAD FIELD + PSEUDOINVERSE =====
% Compressed lead field: column k = sensor topography of eigenmode k
L_tilde = Gain * Phi;                      % [nch x K]

% SVD-based pseudoinverse (correct for K<=nch and K>=nch; no Gram squaring)
[U, S, V] = svd(L_tilde, 'econ');
s = diag(S);

% Rank-safe singular-value floor (MATLAB pinv default if Tol not supplied)
if isempty(Tol)
    if isempty(s)
        Tol = 0;
    else
        Tol = max(size(L_tilde)) * eps(max(s));
    end
end
isKeep = (s > Tol);
sinv = zeros(size(s));
sinv(isKeep) = 1 ./ s(isKeep);

Kernel = V * diag(sinv) * U';              % [K x nch]

%% ===== DIAGNOSTICS =====
rankEff = sum(isKeep);
Info = struct();
Info.CompressedLF   = L_tilde;
Info.SingularValues = s;
Info.Rank           = rankEff;
if rankEff >= 1
    Info.ConditionNumber = s(1) / s(rankEff);
else
    Info.ConditionNumber = Inf;
end
Info.Tol    = Tol;
Info.nModes = size(Phi, 2);
end
