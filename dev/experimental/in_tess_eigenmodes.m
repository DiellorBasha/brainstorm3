function [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile)
% IN_TESS_EIGENMODES: Load precomputed Laplace-Beltrami eigenmodes from a surface file.
%
% USAGE:  [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile)
%
% DESCRIPTION:
%     Loads the embedded Eigenmodes field from a Brainstorm surface file via
%     the canonical loader in_tess_bst. Returns [] and false if eigenmodes
%     have not been computed for this surface.
%
% SEE ALSO: out_tess_eigenmodes, tess_eigenmodes, in_tess_bst

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

Eigenmodes = [];
isComputed = false;

% Load via the canonical Brainstorm loader. isComputeMissing=0 so a
% frequently-called read does not recompute curvature/normals (those recomputes
% are gated on isComputeMissing at in_tess_bst.m:128,135,142,149).
TessMat = in_tess_bst(SurfaceFile, 0);

if ~isfield(TessMat, 'Eigenmodes') || isempty(TessMat.Eigenmodes)
    return;
end

Eigenmodes = TessMat.Eigenmodes;
if isfield(Eigenmodes, 'Vectors') && isa(Eigenmodes.Vectors, 'single')
    Eigenmodes.Vectors = double(Eigenmodes.Vectors);
end
% Backward compatibility: older files have no per-hemisphere metadata.
% Treat them as a single component (rank = column index).
nK = size(Eigenmodes.Vectors, 2);
if ~isfield(Eigenmodes, 'Component') || isempty(Eigenmodes.Component)
    Eigenmodes.Component = ones(nK, 1);
end
if ~isfield(Eigenmodes, 'CompRank') || isempty(Eigenmodes.CompRank)
    Eigenmodes.CompRank = (1:nK)';
end
% Backfill the canonical global eigenvalue order for legacy files (pre-Order).
if ~isfield(Eigenmodes, 'Order') || isempty(Eigenmodes.Order) || numel(Eigenmodes.Order) ~= nK
    [~, Eigenmodes.Order] = sort(double(Eigenmodes.Values(:)), 'ascend');
end
if ~isfield(Eigenmodes, 'nComponents') || isempty(Eigenmodes.nComponents)
    Eigenmodes.nComponents = max(Eigenmodes.Component(:));
end
isComputed = true;
end
