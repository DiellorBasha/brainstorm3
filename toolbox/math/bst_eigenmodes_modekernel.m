function ModeKernel = bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, ImagingKernel, nModes)
% BST_EIGENMODES_MODEKERNEL: Build the mode kernel mapping signals -> mode coefficients.
%
% USAGE:  ModeKernel = bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, ImagingKernel, nModes)
%         Projector  = bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, [], nModes)
%
% DESCRIPTION:
%     Folds the M-weighted eigenmode projection into the imaging kernel, so that
%     applying the result to sensor data yields mode coefficients directly:
%
%         ModeKernel = Phi(:,1:nModes)' * M * ImagingKernel       % [nModes x nChannels]
%         coeff(t)   = ModeKernel * F(t)                          % = Phi'*M*(K*F(t))
%
%     With an empty ImagingKernel, returns the bare projector Phi(:,1:nModes)' * M
%     ([nModes x nVertices]) to apply to a full source matrix (ImageGridAmp).
%
% INPUTS:
%     Eigenmodes    : struct with field .Vectors [nVertices x nModesAll]
%     MassMatrix    : [nVertices x nVertices] sparse mass matrix
%     ImagingKernel : [nVertices x nChannels] imaging kernel (same orientation as
%                     Results.ImagingKernel), or [] for full sources (ImageGridAmp)
%     nModes        : number of leading modes to keep (clamped to available count)
%
% SEE ALSO: bst_eigenmodes_project, tess_laplacian, in_tess_eigenmodes

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

%% ===== PARSE INPUTS =====
Phi = double(Eigenmodes.Vectors);     % [nV x nModesAll]
nV        = size(Phi, 1);
nModesAll = size(Phi, 2);
if (nargin < 4) || isempty(nModes)
    nModes = nModesAll;
end
nModes = max(1, min(nModes, nModesAll));
Phi = Phi(:, 1:nModes);               % [nV x nModes]

%% ===== VALIDATE =====
if (size(MassMatrix,1) ~= nV) || (size(MassMatrix,2) ~= nV)
    error('MassMatrix must be %dx%d.', nV, nV);
end
if (nargin >= 3) && ~isempty(ImagingKernel) && (size(ImagingKernel,1) ~= nV)
    error('ImagingKernel must have %d rows (one per vertex).', nV);
end

%% ===== BUILD MODE KERNEL =====
% Absorb the mass matrix into the kernel once (left-to-right): unlike
% bst_eigenmodes_project (which keeps M*Data sparse-by-dense per call), here the
% projector Phi'*M is precomputed and reused at every time point.
Projector = Phi' * MassMatrix;        % [nModes x nV]
if (nargin < 3) || isempty(ImagingKernel)
    ModeKernel = Projector;                            % [nModes x nV]
else
    ModeKernel = Projector * double(ImagingKernel);    % [nModes x nChannels]
end
end
