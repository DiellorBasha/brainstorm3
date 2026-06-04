function [ConnEig, isComputed] = in_tess_conn_eigenmodes(SurfaceFile)
% IN_TESS_CONN_EIGENMODES: Load connection-Laplacian eigenmodes from a surface file.
%
% USAGE:  [ConnEig, isComputed] = in_tess_conn_eigenmodes(SurfaceFile)
%
% DESCRIPTION:
%     Loads the embedded ConnEigenmodes field via in_tess_bst. Returns [] and
%     false if connection eigenmodes have not been computed for this surface.
%
% SEE ALSO: out_tess_conn_eigenmodes, tess_conn_eigenmodes, in_tess_bst

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

ConnEig    = [];
isComputed = false;

% isComputeMissing=0: a frequent read must not trigger curvature/normal recompute.
TessMat = in_tess_bst(SurfaceFile, 0);

if ~isfield(TessMat, 'ConnEigenmodes') || isempty(TessMat.ConnEigenmodes)
    return;
end

ConnEig = TessMat.ConnEigenmodes;
if isfield(ConnEig, 'Vectors') && isa(ConnEig.Vectors, 'single')
    ConnEig.Vectors = double(ConnEig.Vectors);   % preserves complex
end
nK = size(ConnEig.Vectors, 2);
if ~isfield(ConnEig, 'Component') || isempty(ConnEig.Component)
    ConnEig.Component = ones(nK, 1);
end
if ~isfield(ConnEig, 'CompRank') || isempty(ConnEig.CompRank)
    ConnEig.CompRank = (1:nK)';
end
if ~isfield(ConnEig, 'Order') || isempty(ConnEig.Order) || numel(ConnEig.Order) ~= nK
    [~, ConnEig.Order] = sort(double(ConnEig.Values(:)), 'ascend');
end
if ~isfield(ConnEig, 'nComponents') || isempty(ConnEig.nComponents)
    ConnEig.nComponents = max(ConnEig.Component(:));
end
isComputed = true;
end
