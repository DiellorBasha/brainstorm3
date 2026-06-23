function [L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel, Param, varargin)
% BST_FACE_LEADFIELD  Face-based UNCONSTRAINED MEG forward operator (single-sphere/Sarvas).
%
% USAGE:
%   [L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel, Param)
%   [L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel, Param, 'BlockSize', 500)
%
% DESCRIPTION:
%   Raw Sarvas MEG gain evaluated at the cortical FACE centroids: three Cartesian
%   (x,y,z) columns per face, [nCh x 3F]. This is a pure FORWARD operator -- the
%   full-unconstrained gain, with no orientation applied. It is the face-domain
%   analogue of the standard vertex leadfield and conforms to the Brainstorm
%   convention: the forward stores the unconstrained [nCh x 3N] gain plus the per-
%   source orientation (FaceGeom.Normals -> GridOrient); any orientation constraint
%   (constrained / loose / free) is an inverse-side choice applied DOWNSTREAM via
%   bst_gain_orient(Gain, GridOrient), exactly like the vertex pipeline.
%
%   Geometry is sourced entirely from the canonical manifold_ node (tess_manifold):
%     FaceGeom.Centroids  face centroid positions       (Embedded.face.centroid)
%     FaceGeom.Normals    unit face normals, OUTWARD     (Embedded.face.normal, oriented)
%     FaceGeom.Areas      face areas [m^2]               (Embedded.face.area)
%   The manifold normals are consistently wound but gauge-signed; they are oriented
%   OUTWARD here via the divergence theorem (one global sign per hemisphere), so
%   FaceGeom.Normals follows the physical/outward convention. Faces are returned in
%   the manifold GlobalFaces order (the order the eigen/operator nodes expect).
%
% INPUTS:
%   SurfaceFile  Brainstorm cortex surface file (parent of the manifold_ node)
%   Channel      Brainstorm channel structure (MEG channels)
%   Param        Sphere parameters, one per channel (from bst_headmodeler)
%
% OPTIONS (name-value):
%   'BlockSize'  number of faces per Sarvas call (default 500)
%   'Mode'       'unconstrained' (default and only behavior). 'constrained'/'loose'
%                are NOT a forward concern -- apply bst_gain_orient on the returned
%                gain with FaceGeom.Normals instead.
%
% OUTPUTS:
%   L_face   [nCh x 3*nF]  raw Sarvas gain, 3 Cartesian columns per face
%   FaceGeom struct: .Centroids [nF x 3], .Normals [nF x 3] (outward), .Areas [nF x 1]
%
% SEE ALSO: bst_meg_sph, bst_gain_orient, tess_manifold, bst_dirac
%           dev/references/face_based_source_model.md
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

%% ── Parse options ────────────────────────────────────────────────────────
BlockSize = 500;
Mode      = 'unconstrained';
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'blocksize', BlockSize = varargin{k+1};
        case 'mode',      Mode      = lower(varargin{k+1});
    end
end
if ~strcmpi(Mode, 'unconstrained')
    error(['bst_face_leadfield only produces the UNCONSTRAINED forward gain. ' ...
           '''%s'' is an orientation constraint (an inverse-side concern): apply ' ...
           'bst_gain_orient(L_face, FaceGeom.Normals) on the returned gain instead.'], Mode);
end

%% ── Face geometry from the canonical manifold backbone ───────────────────
M = tess_manifold(SurfaceFile);                         % find-or-load-or-create
TessMat = in_tess_bst(SurfaceFile);
nF = size(TessMat.Faces, 1);
x_f   = zeros(nF, 3);  n_hat = zeros(nF, 3);  A_f = zeros(nF, 1);
for hh = 1:numel(M.Embedded)
    E  = M.Embedded(hh);
    gf = double(E.GlobalFaces);
    x_f(gf,:)   = E.face.centroid;                      % barycentric centroid (exact)
    A_f(gf)     = E.face.area;                          % canonical face area (schemaVersion>=2)
    n_hat(gf,:) = i_orient_outward(E.face.normal, E.face.centroid, E.face.area);  % outward
end

FaceGeom = struct('Centroids', x_f, 'Normals', n_hat, 'Areas', A_f);

%% ── Raw Sarvas gain at face centroids: [nCh x 3F], no projection ─────────
nCh    = numel(Channel);
L_face = zeros(nCh, 3*nF);
for ib = 1:ceil(nF / BlockSize)
    iF   = ((ib-1)*BlockSize + 1) : min(ib*BlockSize, nF);
    cols = reshape([(3*iF-2); (3*iF-1); (3*iF)], 1, []);
    L_face(:, cols) = bst_meg_sph(x_f(iF,:)', Channel, Param);   % raw [nCh x 3*nB]
end
end

%% ────────────────────────────────────────────────────────────────────────
function Nf = i_orient_outward(Nf, Cf, Af)
% Globally flip the consistently-wound manifold face normals to point OUTWARD. For a
% closed surface the divergence theorem gives sum_f Af * n_f . (c_f - c0) = 3*Volume,
% positive iff n_f is outward. One robust global sign per (per-hemisphere) closed submesh.
    c0 = mean(Cf, 1);
    if sum(Af .* sum(Nf .* (Cf - c0), 2)) < 0
        Nf = -Nf;
    end
end
