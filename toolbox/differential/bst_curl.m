function curlField = bst_curl(V, ManifoldMat)
% BST_CURL: Curl (scalar vorticity) of a tangent vector field on the cortical manifold.
%
% On a 2-D surface the de Rham complex is L0 -d0-> L1 -d1-> L2 (scalar -> tangent vector ->
% scalar). Unlike R^3 there is no vector-valued curl: the curl of a tangent field is the SCALAR
% perpendicular vorticity. On the surface curl and divergence are the same operator up to the
% surface Hodge star -- the 90-degree rotation in the tangent plane (the "i"). So the vorticity is
% the divergence of the rotated field:
%   curl V = -div( N x V ) ,   a per-VERTEX scalar,
% where N is the outward face normal (N x V rotates each tangent face vector by 90 degrees). This
% reuses the STABLE adjoint divergence (bst_divergence) and so avoids *1^-1 (the metric flat is
% unstable on obtuse triangles). Because rot(grad f) is divergence-free, curl(grad f) ~ 0. The
% sign convention is CCW-positive seen from OUTSIDE (outward N). v1 handles a TANGENT face field
% (the 3D / normal-bearing case stays with bst_helmholtz). I/O-free: caller passes the loaded node.
%
% USAGE:
%   curlField = bst_curl(V, ManifoldMat)
%
% INPUTS:
%   - V           : tangent vector field, NATIVE per-FACE ambient [3nF x nT] (interleaved x,y,z
%                   per face; e.g. the output of bst_gradient)
%   - ManifoldMat : a manifold_ node (tess_manifold(Surf)) whose 1x2 DEC group + Embedded facet
%                   carry the per-hemisphere operators and face normals.
%
% OUTPUTS:
%   - curlField   : per-VERTEX scalar vorticity field [nV x nT]
%
% SEE ALSO: bst_operators, bst_gradient, bst_divergence, tess_manifold (DEC group)

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

    if ~isfield(ManifoldMat, 'DEC') || isempty(ManifoldMat.DEC) || ~isfield(ManifoldMat.DEC, 'sharp')
        error('bst_curl:noDEC', 'The manifold node has no DEC operators (recompute tess_manifold).');
    end
    DEC = ManifoldMat.DEC;
    Emb = ManifoldMat.Embedded;
    Vrot = zeros(size(V));
    for hh = 1:numel(DEC)
        if isempty(DEC(hh).sharp), continue; end
        fH = double(DEC(hh).GlobalFaces(:));
        % outward face normal (manifold normals are gauge-signed -> orient via divergence theorem)
        Nf = Emb(hh).face.normal;  Cf = Emb(hh).face.centroid;  Af = Emb(hh).face.area;
        c0 = mean(Cf, 1);
        if sum(Af .* sum(Nf .* (Cf - c0), 2)) < 0, Nf = -Nf; end
        % rotate each tangent face vector 90 deg in-plane: Vrot = N x V
        vx = V(3*fH-2, :);  vy = V(3*fH-1, :);  vz = V(3*fH, :);
        nx = Nf(:,1);  ny = Nf(:,2);  nz = Nf(:,3);
        Vrot(3*fH-2, :) = ny.*vz - nz.*vy;
        Vrot(3*fH-1, :) = nz.*vx - nx.*vz;
        Vrot(3*fH,   :) = nx.*vy - ny.*vx;
    end
    curlField = -bst_divergence(Vrot, ManifoldMat);     % curl V = -div(N x V)
end
