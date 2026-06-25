function Wfv = bst_face2vertex(Faces, FaceArea)
% BST_FACE2VERTEX: Area-weighted face->vertex averaging map Wfv [nVh x nFh].
% Wfv * x_face gives, per vertex, the area-weighted mean of its incident faces' values.
% Ported from the bst_helmholtz flat-covariant decomposition; shared by bst_divergence,
% bst_curl, and process_helmholtz.
%
% USAGE:  Wfv = bst_face2vertex(Faces, FaceArea)
%   Faces    : [nFh x 3] LOCAL vertex indices (1..nVh) per face
%   FaceArea : [nFh x 1] face areas
%
% Authors: Diellor Basha, 2026

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

    nFh = size(Faces, 1);  nVh = max(Faces(:));
    I3 = [Faces(:,1); Faces(:,2); Faces(:,3)];  J3 = [(1:nFh)'; (1:nFh)'; (1:nFh)'];
    Wfv = sparse(I3, J3, repmat(FaceArea, 3, 1), nVh, nFh);
    Wfv = spdiags(1 ./ max(sum(Wfv, 2), eps), 0, nVh, nVh) * Wfv;
end
