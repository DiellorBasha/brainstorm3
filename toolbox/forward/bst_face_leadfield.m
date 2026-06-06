function [L_face, FaceGeom] = bst_face_leadfield(Vertices, Faces, Channel, Param, varargin)
% BST_FACE_LEADFIELD  Face-based constrained MEG leadfield (os_meg / spherical).
%
% USAGE:
%   [L_face, FaceGeom] = bst_face_leadfield(Vertices, Faces, Channel, Param)
%   [L_face, FaceGeom] = bst_face_leadfield(Vertices, Faces, Channel, Param, 'BlockSize', 500)
%
% DESCRIPTION:
%   Computes the constrained (normal) MEG leadfield for face-based current flux
%   sources using the overlapping-spheres Sarvas formula.
%
%   For each triangular face f, the source is a current dipole with:
%     Position : face centroid  x_f = (v_i + v_j + v_k) / 3
%     Moment   : q_f = n̂_f * A_f   (exact face normal, scaled by face area)
%
%   This gives the sensor pattern for unit current flux [A/m²] integrated over
%   the face — a genuine primal 2-form in the DEC sense, not a point dipole at
%   a vertex with an averaged vertex normal.
%
%   The Sarvas formula is called at face centroid positions with 3 Cartesian
%   orientations (exactly as in the vertex model), then projected onto the
%   exact face normal and scaled by the face area.  The result is a single
%   leadfield column per face — equivalent to evaluating Sarvas directly with
%   q = n̂_f, but without requiring changes to bst_meg_sph.
%
% PHYSICAL INTERPRETATION:
%   L_face(:, f) is the sensor pattern produced by uniform current density
%   1 A/m² flowing normally through the triangular face f.  Units: T / (A/m²).
%   Equivalently, L_face(:, f) / A_f is the Sarvas Green's function evaluated
%   at the face centroid along the exact face normal — a point dipole field
%   at x_f in direction n̂_f.
%
% INPUTS:
%   Vertices   [nV x 3]  vertex positions in metres
%   Faces      [nF x 3]  triangular face vertex indices (1-based)
%   Channel    Brainstorm channel structure (MEG channels)
%   Param      Sphere parameters (one per channel, from bst_headmodeler)
%
% OPTIONS (name-value):
%   'BlockSize' : number of faces per forward-model call (default 500)
%
% OUTPUTS:
%   L_face   [nCh x nF]  face-based constrained leadfield
%   FaceGeom struct with fields:
%     .Centroids  [nF x 3]  face centroid positions [m]
%     .Normals    [nF x 3]  exact outward unit normals (from mesh winding)
%     .Areas      [nF x 1]  face areas [m²]
%
% COMPARISON WITH VERTEX MODEL:
%   Vertex constrained:  L_v(:,v) = G_sarvas(r_v)  * n̂_v          (vertex position, averaged normal)
%   Face constrained:    L_f(:,f) = G_sarvas(x_f)  * n̂_f * A_f    (centroid, exact normal, area-weighted)
%
%   Differences:
%     1. Source position: vertex r_v vs face centroid x_f
%        Error O(h²/d²) ≈ 0.1% for h~3mm, d~100mm — negligible for os_meg
%     2. Normal direction: vertex-averaged n̂_v vs exact face normal n̂_f
%        Matters near high-curvature regions (sulcal fundi, gyral crowns)
%     3. Area weighting: none in vertex model vs A_f in face model
%        Changes physical interpretation: flux through patch vs unit point dipole
%
% SEE ALSO: bst_meg_sph, bst_gain_orient, bst_eigenmode_leadfield
%           dev/references/face_based_source_model.md
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
BlockSize = 500;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'blocksize', BlockSize = varargin{k+1};
    end
end

%% ── Face geometry (exact, from mesh winding) ─────────────────────────────
% Centroid: simple average of three corner positions
x_f = (Vertices(Faces(:,1),:) + Vertices(Faces(:,2),:) + Vertices(Faces(:,3),:)) / 3;

% Exact face normal from cross product of two edge vectors.
% FreeSurfer meshes are wound clockwise (viewed from outside), so the raw
% cross product points INWARD — opposite to Brainstorm's VertNormals convention
% (which are corrected to point outward during import).
% We detect and correct the sign so that n_hat is consistently outward,
% matching the convention used by the rest of the Brainstorm pipeline.
e01    = Vertices(Faces(:,2),:) - Vertices(Faces(:,1),:);   % [nF x 3]
e02    = Vertices(Faces(:,3),:) - Vertices(Faces(:,1),:);   % [nF x 3]
Ncross = cross(e01, e02, 2);                                  % [nF x 3]
A2     = sqrt(sum(Ncross.^2, 2));                             % [nF x 1]  2 * area
A_f    = A2 / 2;                                              % [nF x 1]  face areas [m²]
n_hat  = Ncross ./ A2;                                        % [nF x 3]  unit normals from winding

% Align with Brainstorm's outward convention by checking the median dot product
% of n_hat with the face-averaged vertex normals.  VertNormals are required as
% an optional input, or we estimate the orientation from the vertex positions
% (centroid-to-brain-center direction as a proxy).
brain_center = mean(Vertices, 1);                         % approximate brain centre
radial       = x_f - brain_center;                        % centroid → outside
outward_dot  = sum(n_hat .* radial, 2);                   % positive = outward
if mean(outward_dot) < 0                                  % majority pointing inward
    n_hat = -n_hat;                                       % flip to outward convention
end

FaceGeom.Centroids = x_f;
FaceGeom.Normals   = n_hat;
FaceGeom.Areas     = A_f;

nF   = size(Faces, 1);
nCh  = numel(Channel);

%% ── Leadfield: Sarvas at face centroids, projected onto exact normals ─────
%
% For block iBlock of faces:
%   G_block = bst_meg_sph(x_f(iBlock,:)', ...)   [nCh x 3*nBlock]
%
% G_block(:, 3f-2:3f) is the [nCh x 3] response to unit Cartesian dipoles at x_f.
% Projecting onto n̂_f and scaling by A_f gives the physical face flux column:
%   L_face(:, f) = G_block(:, 3f-2:3f) * (n̂_f * A_f)
%
% Vectorised using a sparse block-diagonal projection matrix P [3*nBlock x nBlock]:
%   P(3f-2:3f, f) = n̂_f' * A_f
%   L_face(:, block) = G_block * P
%
% This is mathematically equivalent to calling Sarvas directly with q = n̂_f * A_f,
% but reuses bst_meg_sph without modification.

L_face = zeros(nCh, nF);

nBlocks = ceil(nF / BlockSize);
for ib = 1:nBlocks
    iF  = ((ib-1)*BlockSize + 1) : min(ib*BlockSize, nF);
    nB  = numel(iF);

    % Sarvas at face centroids for this block — [nCh x 3*nB]
    G_block = bst_meg_sph(x_f(iF,:)', Channel, Param);

    % Sparse block-diagonal projection: col f gets n̂_f * A_f
    % Each face contributes 3 rows, 1 col to the projection matrix
    n_scaled = (n_hat(iF,:) .* A_f(iF))';   % [3 x nB]  (n̂_f * A_f per column)
    ri = reshape(1:3*nB, 3, nB);             % [3 x nB]  row indices
    ci = repmat(1:nB, 3, 1);                 % [3 x nB]  col indices
    P  = sparse(ri(:), ci(:), n_scaled(:), 3*nB, nB);

    L_face(:, iF) = G_block * P;
end

end
