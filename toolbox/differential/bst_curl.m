function curlField = bst_curl(V, ManifoldMat, varargin)
% BST_CURL: Curl (scalar vorticity) of a tangent or ambient vector field on the cortical manifold.
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
% sign convention is CCW-positive seen from OUTSIDE (outward N).
%
% Two field types are supported:
%   TANGENT: v1 per-FACE field [3nF x nT] (interleaved x,y,z per face; e.g. output of
%            bst_gradient). Uses -div(N x V) via the stable DEC adjoint.
%            Call: bst_curl(V, ManifoldMat).
%   AMBIENT: 3D R^3 per-VERTEX field [3nV x nT] (interleaved x,y,z per vertex; e.g. a Dirac
%            source vector). Returns the Hodge vorticity (Dirac w-part) from bst_helmholtz.
%            No mean-curvature coupling for curl (vorticity is purely intrinsic).
%            Call: bst_curl(V, ManifoldMat, 'Ambient', Surf, Dir, LBO).
% I/O-free: caller passes the loaded node (and operator nodes for ambient branch).
%
% USAGE:
%   curlField = bst_curl(V, ManifoldMat)
%   curlField = bst_curl(V, ManifoldMat, 'Ambient', Surf, Dir, LBO)
%
% INPUTS:
%   - V           : [TANGENT] per-FACE ambient [3nF x nT]; [AMBIENT] per-VERTEX [3nV x nT]
%   - ManifoldMat : a manifold_ node (tess_manifold(Surf)) whose 1x2 DEC group + Embedded facet
%                   carry the per-hemisphere operators and face normals.
%   - 'Ambient'   : (optional) flag to dispatch to the ambient branch
%   - Surf        : (ambient only) loaded tessellation struct (in_tess_bst output)
%   - Dir         : (ambient only) Dirac operator node (bst_get_operator_node(Surf,'Dirac'))
%   - LBO         : (ambient only) LBO operator node (bst_get_operator_node(Surf,'Laplace-Beltrami'))
%
% OUTPUTS:
%   - curlField   : per-VERTEX scalar vorticity field [nV x nT]
%
% SEE ALSO: bst_operators, bst_gradient, bst_divergence, bst_helmholtz, tess_manifold (DEC group)

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

    % ----- ambient (3nV) branch: Hodge vorticity (Dirac w-part) -----
    if ~isempty(varargin) && strcmpi(varargin{1}, 'Ambient')
        Surf = varargin{2};  Dir = varargin{3};  LBO = varargin{4};
        H = bst_helmholtz('Decompose', {Dir, LBO}, ManifoldMat, Surf, V);
        curlField = H.Curl;
        return;
    end
    % ----- tangent (3nF) branch: existing -div(N x V) (UNCHANGED) -----
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
