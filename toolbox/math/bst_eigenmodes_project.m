function [Coeffs, Reconstructed] = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix, varargin)
% BST_EIGENMODES_PROJECT: Project scalar data onto a precomputed eigenmode basis.
%
% USAGE:  Coeffs = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix)
%         Coeffs = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix, 'ModeRange', [1 100])
%         [Coeffs, Reconstructed] = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix, ...)
%
% INPUTS:
%     Eigenmodes : struct from in_tess_eigenmodes (uses .Vectors)
%     Data       : [nVertices x nTime] scalar field(s) on the mesh
%     MassMatrix : [nVertices x nVertices] sparse mass matrix (from tess_laplacian)
%
% OPTIONS:
%     ModeRange  : [1 x 2] mode index range for reconstruction (default: all)
%
% SEE ALSO: tess_eigenmodes, in_tess_eigenmodes, bst_eigenmodes_filter, tess_laplacian

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
ModeRange = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'moderange', ModeRange = varargin{i+1};
    end
end

Phi    = double(Eigenmodes.Vectors);   % [nV x nModes]
nV     = size(Phi, 1);
nModes = size(Phi, 2);

%% ===== VALIDATE =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but eigenmodes have %d vertices.', size(Data, 1), nV);
end
if (size(MassMatrix, 1) ~= nV) || (size(MassMatrix, 2) ~= nV)
    error('MassMatrix must be %dx%d.', nV, nV);
end

%% ===== FORWARD TRANSFORM (all modes, stored order) =====
Coeffs = manifold_ft(Phi, MassMatrix, Data);   % [nModes x nTime]

%% ===== RECONSTRUCT over the CANONICAL mode range (if requested) =====
if nargout >= 2
    if isempty(ModeRange); ModeRange = [1, nModes]; end
    ModeRange(1) = max(1, ModeRange(1));
    ModeRange(2) = min(nModes, ModeRange(2));
    if ModeRange(2) < ModeRange(1)
        error('Empty mode range [%d, %d] (have %d modes).', ModeRange(1), ModeRange(2), nModes);
    end
    if isfield(Eigenmodes,'Order') && ~isempty(Eigenmodes.Order)
        Order = double(Eigenmodes.Order(:));
    else
        [~, Order] = sort(double(Eigenmodes.Values(:)), 'ascend');
    end
    iSel = Order(ModeRange(1):ModeRange(2));
    Reconstructed = manifold_ift(Phi(:, iSel), Coeffs(iSel, :));   % [nV x nTime]
end
end
