function M = bst_eigenmodes_harmonic(L, Phi, iW)
% BST_EIGENMODES_HARMONIC: Unregularized whitened eigenmode imaging kernel.
%
% USAGE:  M = bst_eigenmodes_harmonic(L, Phi, iW)
%
% DESCRIPTION:
%     Returns the eigenmode-space "Harmonic" kernel
%         M = pinv(iW * L * Phi) * iW            [K x nCh]
%     i.e. the rank-safe (SVD) pseudoinverse of the whitened compressed lead
%     field, applied to the whitener so it maps RAW recordings to eigenmode
%     coefficients. There is no source/inverse regularization (no Tikhonov, no
%     prior) -- only the rank-safe singular-value floor of bst_eigenmodes_transform.
%     The vertex-space source map is Phi*M; the eigenmode time series is M*Data.
%
% INPUTS:
%     L   : [nCh x nVert] constrained (fixed-orientation) lead field, good channels.
%     Phi : [nVert x K]   eigenmode matrix (caller truncates to K modes).
%     iW  : [nCh x nCh]   noise whitener (pass eye(nCh) for no whitening).
%
% OUTPUT:
%     M   : [K x nCh] harmonic eigenmode kernel.
%
% SEE ALSO: bst_eigenmodes_transform, bst_inverse_eigenmodes

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

% Rank-safe pseudoinverse of the whitened compressed lead field; folding in iW
% lets M pre-whiten raw recordings on the fly (M*d = pinv(iW*L*Phi)*(iW*d)).
[Kt, ~] = bst_eigenmodes_transform(iW * L, Phi);   % pinv(iW*L*Phi)  [K x nCh]
M = Kt * iW;                                        % [K x nCh]
end
