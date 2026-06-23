function [L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel, Param, varargin)
% BST_FACE_LEADFIELD  Face-based MEG leadfield using the trivial-connection frame.
%
% USAGE:
%   [L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel, Param)
%   [L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel, Param, 'BlockSize', 500)
%   [L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel, Param, 'Mode', 'loose')
%
% DESCRIPTION:
%   Face-based constrained (or loose) MEG leadfield.  The tangent frame at each
%   face comes from the FreeSurfer trivial-connection solve (tess_tangents):
%
%     {U_f, V_f}  — smooth, globally consistent per-face tangent vectors
%     n̂_f = U_f × V_f  — exact face normal, sign-consistent by winding + orientation
%
%   This means all three frame vectors are available for free and are already
%   in the correct gauge for downstream wave analysis, DEC gradient computation,
%   and the connection-Laplacian eigenmode framework.  No cross-product of edge
%   vectors, no vertex-normal averaging, no sign-correction heuristic.
%
%   Source at face f:
%     Constrained:  q_f = s_n · n̂_f · A_f          (normal only, 1 unknown)
%     Loose:        q_f = s_n · n̂_f · A_f
%                       + s1 · U_f               (+ 2 tangential, 3 unknowns)
%                       + s2 · V_f
%
%   The Sarvas formula is called at face centroid positions, then projected
%   onto each frame direction. Equivalent to calling Sarvas with q = n̂_f
%   (or U_f, V_f) directly — zero approximation in the direction.
%   O(h²/d²) ≈ 0.1% approximation from using centroid vs integrating over face.
%
% MODES:
%   'unconstrained'         : FULL-3D, three Cartesian (x,y,z) columns per face =
%                             the raw Sarvas gain at the centroid [nCh x 3F]. No
%                             frame, no projection, no tess_tangents -- exact parity
%                             with the vertex leadfield. Geometry from tess_manifold.
%                             The rigorous path for the full-unconstrained inverse.
%   'constrained' (default) : one leadfield column per face along n̂_f, scaled A_f
%   'loose'                 : three columns per face [n̂_f*A_f, U_f, V_f]
%                             Normal column scaled by A_f (flux units);
%                             tangential columns unit (direction field units).
%                             (constrained/loose use the legacy tess_tangents frame.)
%
% INPUTS:
%   SurfaceFile  Brainstorm surface file (must have Reg.Sphere for tess_tangents)
%   Channel      Brainstorm channel structure (MEG channels)
%   Param        Sphere parameters (one per channel, from bst_headmodeler)
%
% OPTIONS (name-value):
%   'BlockSize'  number of faces per Sarvas call (default 500)
%   'Mode'       'constrained' | 'loose'  (default 'constrained')
%
% OUTPUTS:
%   L_face   [nCh x nF]    constrained: one column per face
%            [nCh x 3*nF]  loose: three columns per face (n, U, V interleaved)
%   FaceGeom struct:
%     .Centroids  [nF x 3]  face centroid positions [m]
%     .Normals    [nF x 3]  n̂_f = U_f × V_f  (outward, trivial-connection sign)
%     .Areas      [nF x 1]  face areas [m²]
%     .U          [nF x 3]  trivial-connection tangent e1 (smooth, consistent gauge)
%     .V          [nF x 3]  trivial-connection tangent e2 = n̂_f × U_f
%
% RELATIONSHIP TO VERTEX MODEL:
%   Vertex constrained:  L_v(:,v) = G(r_v) · n̂_v
%     — point at vertex, averaged normal, no area, nxr CW sign
%   Face constrained:    L_face(:,f) = G(x_f) · n̂_f · A_f
%     — centroid, exact trivial-connection normal, area-weighted, outward sign
%
% SEE ALSO: tess_tangents, bst_meg_sph, bst_gain_orient, bst_eigenmode_leadfield
%           dev/references/face_based_source_model.md
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
BlockSize = 500;
Mode      = 'constrained';
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'blocksize', BlockSize = varargin{k+1};
        case 'mode',      Mode      = lower(varargin{k+1});
    end
end
doLoose = strcmpi(Mode, 'loose');

%% ── FULL-UNCONSTRAINED mode (the rigorous full-3D path) ──────────────────
% Raw Sarvas gain at face centroids = [nCh x 3F], three Cartesian (x,y,z) columns
% per face. No frame, no projection, no tess_tangents (deprecated) -- exact parity
% with the vertex leadfield. Geometry comes from tess_manifold (centroids + unit
% gauge-consistent normals); normals are GridOrient display only, never projected.
if strcmpi(Mode, 'unconstrained')
    [L_face, FaceGeom] = i_unconstrained(SurfaceFile, Channel, Param, BlockSize);
    return;
end

%% ── Trivial-connection face frame (the only frame we need) ───────────────
% tess_tangents returns per-face vectors from the FreeSurfer trivial-connection
% solve: U_f (e1) and V_f (e2) span the face tangent plane and are smooth +
% globally consistent across the hemisphere.  The face normal is their cross
% product — exact, sign-consistent, no averaging.
[U_f, V_f] = tess_tangents(SurfaceFile, 'NoSave', 1);   % [nF x 3] each

% Face normal from the trivial-connection frame orientation
n_hat = cross(U_f, V_f, 2);                              % [nF x 3]
n_nrm = sqrt(sum(n_hat.^2, 2));
n_hat = n_hat ./ n_nrm;                                  % unit normals

% Face areas and centroids from mesh geometry
TessMat  = in_tess_bst(SurfaceFile);
Vertices = TessMat.Vertices;
Faces    = double(TessMat.Faces);
x_f      = (Vertices(Faces(:,1),:) + Vertices(Faces(:,2),:) + Vertices(Faces(:,3),:)) / 3;

edge_a = Vertices(Faces(:,2),:) - Vertices(Faces(:,1),:);
edge_b = Vertices(Faces(:,3),:) - Vertices(Faces(:,1),:);
A_f    = sqrt(sum(cross(edge_a, edge_b, 2).^2, 2)) / 2;   % [nF x 1]

FaceGeom.Centroids = x_f;
FaceGeom.Normals   = n_hat;
FaceGeom.Areas     = A_f;
FaceGeom.U         = U_f;
FaceGeom.V         = V_f;

nF  = size(Faces, 1);
nCh = numel(Channel);

%% ── Leadfield: Sarvas at face centroids, projected onto trivial-connection frame
%
% For block iB of faces:
%   G_block = bst_meg_sph(x_f(iB,:)', ...)         [nCh x 3*nB]
%
% Constrained — project onto n̂_f, scale by A_f:
%   L_face(:, f) = G_block(:, 3f-2:3f) * (n̂_f * A_f)
%
% Loose — three columns per face:
%   L_face(:, 3f-2) = G_block(:, 3f-2:3f) * n̂_f * A_f    normal, flux units
%   L_face(:, 3f-1) = G_block(:, 3f-2:3f) * U_f           tangential e1
%   L_face(:, 3f  ) = G_block(:, 3f-2:3f) * V_f           tangential e2
%
% Both computed via a sparse block-diagonal projection matrix.

nCols  = 1 + 2*doLoose;          % 1 for constrained, 3 for loose
L_face = zeros(nCh, nCols * nF);

nBlocks = ceil(nF / BlockSize);
for ib = 1:nBlocks
    iF = ((ib-1)*BlockSize + 1) : min(ib*BlockSize, nF);
    nB = numel(iF);

    G_block = bst_meg_sph(x_f(iF,:)', Channel, Param);   % [nCh x 3*nB]

    if ~doLoose
        % Constrained: one column per face = G * (n̂_f * A_f)
        n_scaled = (n_hat(iF,:) .* A_f(iF))';            % [3 x nB]
        ri = reshape(1:3*nB, 3, nB);
        ci = repmat(1:nB, 3, 1);
        P  = sparse(ri(:), ci(:), n_scaled(:), 3*nB, nB);
        L_face(:, iF) = G_block * P;
    else
        % Loose: three columns per face — interleaved [n̂*A, U, V] per face
        % Column order in output: face 1 → cols 1,2,3; face 2 → cols 4,5,6; ...
        % Projection matrix [3*nB x 3*nB]:
        %   For face f (local index k): rows 3k-2:3k, cols 3k-2 = n̂_f*A_f
        %                                                col  3k-1 = U_f
        %                                                col  3k   = V_f
        proj = zeros(3*nB, 3*nB);
        for k = 1:nB
            rows = (3*k-2):(3*k);
            proj(rows, 3*k-2) = n_hat(iF(k),:)' * A_f(iF(k));   % normal, flux
            proj(rows, 3*k-1) = U_f(iF(k),:)';                    % tangential e1
            proj(rows, 3*k  ) = V_f(iF(k),:)';                    % tangential e2
        end
        iCols = reshape([(3*(iF-1)+1); (3*(iF-1)+2); (3*iF)], 1, []);
        L_face(:, iCols) = G_block * sparse(proj);
    end
end

end


%% ════════════════════════════════════════════════════════════════════════
function [L_face, FaceGeom] = i_unconstrained(SurfaceFile, Channel, Param, BlockSize)
% Full-3D unconstrained face leadfield: raw Sarvas gain at face centroids,
% [nCh x 3F], three Cartesian columns per face. Geometry from tess_manifold.

    % --- face geometry from the canonical manifold backbone (not tess_tangents) ---
    M = tess_manifold(SurfaceFile);                         % find-or-load-or-create
    TessMat  = in_tess_bst(SurfaceFile);
    Faces = double(TessMat.Faces);
    nF = size(Faces, 1);
    x_f   = zeros(nF, 3);  n_hat = zeros(nF, 3);  A_f = zeros(nF, 1);
    for hh = 1:numel(M.Embedded)
        gf = M.Embedded(hh).GlobalFaces;
        x_f(gf,:)   = M.Embedded(hh).face.centroid;         % == barycentric centroid (exact)
        n_hat(gf,:) = M.Embedded(hh).face.normal;           % unit, gauge-consistent orientation
        A_f(gf)     = M.Embedded(hh).face.area;             % canonical face area (schemaVersion>=2)
    end

    FaceGeom = struct('Centroids', x_f, 'Normals', n_hat, 'Areas', A_f);

    % --- raw Sarvas gain: [nCh x 3F], 3 Cartesian (x,y,z) columns per face, NO projection ---
    nCh    = numel(Channel);
    L_face = zeros(nCh, 3*nF);
    for ib = 1:ceil(nF / BlockSize)
        iF = ((ib-1)*BlockSize + 1) : min(ib*BlockSize, nF);
        cols = reshape([(3*iF-2); (3*iF-1); (3*iF)], 1, []);
        L_face(:, cols) = bst_meg_sph(x_f(iF,:)', Channel, Param);   % raw [nCh x 3*nB]
    end
end
