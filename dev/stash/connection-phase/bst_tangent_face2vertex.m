function [Uv, Vv] = bst_tangent_face2vertex(Faces, Uf, VertNormals)
% BST_TANGENT_FACE2VERTEX: Per-face tangent direction -> per-vertex orthonormal frame.
%
% USAGE:  [Uv, Vv] = bst_tangent_face2vertex(Faces, Uf, VertNormals)
%
% DESCRIPTION:
%     Transfers a per-FACE tangent direction Uf (e.g. the trivial-connection e1
%     from tess_tangents) to a per-VERTEX orthonormal tangent frame. Each vertex
%     averages the e1 vectors of its incident faces, projects the result into the
%     vertex tangent plane (removing the vertex-normal component) and renormalizes;
%     Vv = n x Uv completes a right-handed frame. The trivial-connection field is
%     smooth, so plain averaging is valid away from the (few) singularities, where
%     the incident directions cancel and Uv is left as an arbitrary in-plane unit
%     vector (those vertices are the field singularities).
%
% INPUTS:
%     Faces       : [nF x 3] 1-based triangle indices.
%     Uf          : [nF x 3] per-face tangent direction (unit).
%     VertNormals : [nV x 3] per-vertex normals.
%
% OUTPUTS:
%     Uv : [nV x 3] per-vertex first tangent (unit, in the tangent plane).
%     Vv : [nV x 3] per-vertex second tangent = n x Uv.
%
% SEE ALSO: tess_tangents, bst_conn_phase

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

Faces = double(Faces);
nV = size(VertNormals, 1);

% Accumulate incident-face e1 at each vertex (uniform average over incident faces).
ii   = [Faces(:,1); Faces(:,2); Faces(:,3)];   % [3nF x 1]
Urep = [Uf;          Uf;          Uf];          % [3nF x 3]
Usum = zeros(nV, 3);
for d = 1:3
    Usum(:, d) = accumarray(ii, Urep(:, d), [nV, 1]);
end

% Unit vertex normals.
N = VertNormals ./ max(sqrt(sum(VertNormals.^2, 2)), eps);

% Project the averaged direction into the tangent plane and renormalize.
Uv  = Usum - sum(Usum .* N, 2) .* N;
nrm = sqrt(sum(Uv.^2, 2));
% Guard singularities (cancellation): fall back to an arbitrary in-plane axis.
bad = nrm < 1e-9;
if any(bad)
    ref = repmat([1 0 0], sum(bad), 1);
    alt = abs(N(bad,1)) > 0.9;          % avoid degeneracy when n ~ x-axis
    ref(alt, :) = repmat([0 1 0], sum(alt), 1);
    proj = ref - sum(ref .* N(bad,:), 2) .* N(bad,:);
    Uv(bad, :) = proj;
    nrm(bad)   = sqrt(sum(proj.^2, 2));
end
Uv = Uv ./ max(nrm, eps);
Vv = cross(N, Uv, 2);
end
