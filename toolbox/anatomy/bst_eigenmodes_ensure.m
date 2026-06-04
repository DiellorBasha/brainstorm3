function Eig = bst_eigenmodes_ensure(SurfaceFile, nModesPerHemi)
% BST_EIGENMODES_ENSURE: Return a surface's canonical eigenmodes, computing a default
% set if absent. The principled home for the "use the existing eigenmodes of a surface,
% else compute a default" policy.
%
% USAGE:  Eig = bst_eigenmodes_ensure(SurfaceFile)
%         Eig = bst_eigenmodes_ensure(SurfaceFile, nModesPerHemi)   % default 1000
%
% DESCRIPTION:
%     If the surface already carries eigenmodes, they are returned as-is (with the
%     canonical .Order backfilled by in_tess_eigenmodes). Otherwise a default set is
%     computed (nModesPerHemi modes per connected component, barycentric mass, DC mode
%     removed) and stored on the surface. NO repair is attempted: a non-manifold surface
%     raises an error, because repair changes the vertex count and breaks the
%     surface <-> leadfield <-> eigenmode consistency. Remesh to an icosphere instead.
%
% SEE ALSO: tess_eigenmodes, in_tess_eigenmodes, out_tess_eigenmodes, tess_manifold
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

if nargin < 2 || isempty(nModesPerHemi)
    nModesPerHemi = 1000;
end

% Reuse existing eigenmodes when present (canonical .Order backfilled on read).
[Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
if isComputed && ~isempty(Eig)
    return;
end

% Compute a default set. Guard manifoldness first -- no silent repair.
Surf = in_tess_bst(SurfaceFile, 0);
mani = tess_manifold(Surf.Vertices, Surf.Faces);
if isstruct(mani) && isfield(mani, 'isManifold') && ~mani.isManifold
    error('bst_eigenmodes_ensure:NonManifold', ...
        ['Surface %s is non-manifold; eigenmodes require a 2-manifold mesh. ' ...
         'Remesh to an icosphere (or repair manually) and retry.'], SurfaceFile);
end

Eig = tess_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', nModesPerHemi, ...
    'MassType', 'barycentric', 'RemoveDC', 1);
out_tess_eigenmodes(SurfaceFile, Eig);
end
