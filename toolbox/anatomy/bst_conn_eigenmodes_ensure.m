function ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile, nModesPerHemi)
% BST_CONN_EIGENMODES_ENSURE: Return a surface's canonical connection eigenmodes,
% computing a default set if absent. Sibling of bst_eigenmodes_ensure for the
% vector-field (connection-Laplacian) axis.
%
% USAGE:  ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile)
%         ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile, nModesPerHemi)
%
% DESCRIPTION:
%     If the surface already carries ConnEigenmodes, they are returned as-is.
%     Otherwise a default set is computed and stored. When nModesPerHemi is not
%     given, the per-component count is matched to the surface's scalar Eigenmodes
%     axis as round(nModes / nComponents) -- the per-component average, exact for
%     equal-sized components and within one mode otherwise -- ensuring that axis
%     exists first via bst_eigenmodes_ensure. NO repair is attempted: a non-manifold
%     surface raises
%     an error (repair changes the vertex count and breaks surface<->eigenmode
%     consistency). Remesh to an icosphere instead.
%
% SEE ALSO: tess_conn_eigenmodes, in_tess_conn_eigenmodes, out_tess_conn_eigenmodes,
%           bst_eigenmodes_ensure, tess_manifold

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

% Reuse existing connection eigenmodes when present.
[ConnEig, isComputed] = in_tess_conn_eigenmodes(SurfaceFile);
if isComputed && ~isempty(ConnEig)
    return;
end

% Derive the per-component count from the scalar axis if not given (ensuring the
% scalar axis exists first; bst_eigenmodes_ensure reuses it when already present).
if nargin < 2 || isempty(nModesPerHemi)
    sEig = bst_eigenmodes_ensure(SurfaceFile);
    nModesPerHemi = max(1, round(sEig.nModes / max(1, sEig.nComponents)));
end

% Compute a default set. tess_conn_eigenmodes validates the 2-manifold itself
% (no silent repair), raising tess_conn_eigenmodes:NonManifold on a bad mesh.
Surf = in_tess_bst(SurfaceFile, 0);

ConnEig = tess_conn_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', nModesPerHemi);
out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Surf.Vertices, Surf.Faces);
end
