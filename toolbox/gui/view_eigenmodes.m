function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Browse Laplace-Beltrami eigenmodes on a surface.
%
% USAGE:  hFig = view_eigenmodes(SurfaceFile)
%         [Grid, K, Info] = view_eigenmodes('BuildPairedGrid', Eigenmodes)
%
% The viewer displays each component's rank-k mode together (mode k shows both
% hemispheres) as a registered Brainstorm Source result, so the standard colormap
% UI applies; Left/Right arrows step modes.
%
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

if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'BuildPairedGrid')
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: paired display grid (column k = each component's CompRank==k mode) =====
function [Grid, K, Info] = BuildPairedGrid(Eig)
    nV = size(Eig.Vectors, 1);
    nK = size(Eig.Vectors, 2);
    if isfield(Eig, 'CompRank') && ~isempty(Eig.CompRank)
        CompRank  = Eig.CompRank(:);
        Component = Eig.Component(:);
    else
        CompRank  = (1:nK)';
        Component = ones(nK, 1);
    end
    K = max(CompRank);
    Grid = zeros(nV, K);
    for k = 1:K
        cols = find(CompRank == k);
        if ~isempty(cols)
            Grid(:, k) = sum(Eig.Vectors(:, cols), 2);   % disjoint support across components
        end
    end
    Info = struct('K', K, 'Component', Component, 'CompRank', CompRank, 'Values', Eig.Values(:));
end
