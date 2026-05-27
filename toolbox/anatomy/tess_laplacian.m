function [L, M] = tess_laplacian(Vertices, Faces, varargin)
% TESS_LAPLACIAN: Cotangent Laplacian and mass matrix for a triangle mesh.
%
% USAGE:  [L, M] = tess_laplacian(Vertices, Faces)
%         [L, M] = tess_laplacian(Vertices, Faces, 'MassType', 'barycentric')
%
% DESCRIPTION:
%     Assembles the cotangent Laplacian (stiffness matrix) and lumped mass
%     matrix for a triangulated surface, with no external dependencies.
%
%     The cotangent Laplacian L is the standard FEM discretization of the
%     Laplace-Beltrami operator on triangle meshes:
%       L(i,j) = -0.5 * [cot(alpha_ij) + cot(beta_ij)]   for edge (i,j)
%       L(i,i) = -sum_j L(i,j)
%     where alpha_ij and beta_ij are the angles opposite edge (i,j) in the
%     two adjacent triangles. This yields a POSITIVE SEMIDEFINITE matrix
%     (eigenvalues >= 0), with the convention that the Dirichlet energy
%     E(u) = u' * L * u >= 0.
%
%     The mass matrix M defines the L2 inner product on the mesh:
%       <u, v>_M = u' * M * v
%     Required for the generalized eigenvalue problem L*phi = lambda*M*phi.
%
% INPUTS:
%     Vertices : [nVertices x 3] vertex positions
%     Faces    : [nFaces x 3] triangle vertex indices (1-based)
%
% OPTIONS:
%     MassType : Type of mass matrix (default: 'barycentric')
%                'barycentric' — Lumped diagonal, area/3 per vertex per face.
%                                Simple, robust, always positive diagonal.
%                'voronoi'     — Lumped diagonal, Voronoi dual cell areas.
%                                More accurate for irregular meshes but
%                                requires obtuse-triangle handling.
%                'galerkin'    — Consistent FEM (sparse, non-diagonal).
%                                M(i,j) = area/12 for edge (i,j) in face,
%                                M(i,i) = area/6 for vertex i in face.
%     Symmetrize : (logical) Enforce symmetry via (L+L')/2. Default: true.
%     CheckManifold : (logical) Warn if the input mesh is not a clean
%                     2-manifold. Default: false. This function ASSUMES a
%                     clean 2-manifold input (e.g. icosphere-downsampled
%                     cortex); validate upstream with tess_manifold.
%
% OUTPUTS:
%     L : [nVertices x nVertices] sparse positive semidefinite cotangent
%         Laplacian (stiffness matrix)
%     M : [nVertices x nVertices] sparse mass matrix (diagonal for
%         barycentric/voronoi, sparse for galerkin)
%
% MATHEMATICAL REFERENCE:
%     Cotangent weights: For triangle (i,j,k), the edge (i,j) gets weight
%       w_ij = 0.5 * cot(angle at k)
%     where cot(theta) = dot(u,v) / ||cross(u,v)|| with u = Vi-Vk, v = Vj-Vk.
%
%     The full weight for edge (i,j) sums contributions from both adjacent
%     triangles: W_ij = 0.5*[cot(alpha) + cot(beta)].
%
%     For boundary edges (one adjacent triangle), only one cotangent
%     contributes. This implements natural (Neumann) boundary conditions.
%
% SEE ALSO: tess_manifold

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

%% ===== PARSE INPUTS =====
MassType = 'barycentric';
Symmetrize = true;
CheckManifold = false;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'masstype',      MassType = lower(varargin{i+1});
        case 'symmetrize',    Symmetrize = varargin{i+1};
        case 'checkmanifold', CheckManifold = varargin{i+1};
    end
end

% Validate inputs
if (size(Vertices, 2) ~= 3) || (size(Faces, 2) ~= 3)
    error('Vertices and Faces must have 3 columns.');
end
nV = size(Vertices, 1);
Faces = double(Faces);

% Optional manifold check (OFF by default). This operator is recomputed on
% every project/filter call, and tess_eigenmodes is the authoritative gate,
% so a full manifold scan in the inner loop would be wasteful.
if CheckManifold
    [~, ~, isManifold] = tess_manifold(Vertices, Faces, 'Repair', 0, 'Verbose', 0);
    if ~isManifold
        warning('tess_laplacian:NonManifold', ...
            ['Input mesh is not a clean 2-manifold; the cotangent Laplacian may be ' ...
             'inaccurate. Validate with tess_manifold or re-mesh with icosphere downsampling.']);
    end
end

%% ===== COTANGENT LAPLACIAN =====
% For each face (i,j,k), compute cotangent weights for all three edges.
%
% Edge opposite vertex 1 (i.e., edge 2-3):
%   cot(angle at vertex 1) = dot(V2-V1, V3-V1) / ||cross(V2-V1, V3-V1)||
%
% We compute all three angles per face simultaneously.

% Extract vertex positions for each face corner
V1 = Vertices(Faces(:,1), :);  % [nF x 3]
V2 = Vertices(Faces(:,2), :);
V3 = Vertices(Faces(:,3), :);

% Edge vectors
e1 = V2 - V1;  % edge opposite vertex 3... no, edge from V1 to V2
e2 = V3 - V1;  % edge from V1 to V3
e3 = V3 - V2;  % edge from V2 to V3

% Cotangent at vertex 1 (angle between edges e1 and e2):
%   cot = dot(e1, e2) / ||cross(e1, e2)||
% This is the weight for the OPPOSITE edge (2,3).
dot12 = sum(e1 .* e2, 2);
cross12 = cross(e1, e2, 2);
dblArea = sqrt(sum(cross12.^2, 2));  % = 2 * triangle area

% Cotangent at vertex 2 (angle between edges -e1 and e3):
%   Vectors from V2: (V1-V2) = -e1, (V3-V2) = e3
dot_at2 = sum((-e1) .* e3, 2);

% Cotangent at vertex 3 (angle between edges -e2 and -e3):
%   Vectors from V3: (V1-V3) = -e2, (V2-V3) = -e3
dot_at3 = sum((-e2) .* (-e3), 2);

% All three cotangents share the same denominator (||cross|| = 2*Area)
% because the cross product magnitude is the same regardless of which
% vertex the angle is at (it's twice the triangle area).
denom = dblArea;
denom(denom == 0) = eps;  % protect against degenerate triangles

cot1 = dot12 ./ denom;     % cot(angle at V1), weight for edge (V2,V3)
cot2 = dot_at2 ./ denom;   % cot(angle at V2), weight for edge (V1,V3)
cot3 = dot_at3 ./ denom;   % cot(angle at V3), weight for edge (V1,V2)

% Scale by 0.5 (FEM convention)
cot1 = 0.5 * cot1;
cot2 = 0.5 * cot2;
cot3 = 0.5 * cot3;

% Assemble off-diagonal entries of L.
% Each face contributes 3 edges, each edge gets the cotangent of the
% opposite vertex. For an interior edge, both adjacent faces contribute.
%
% L(i,j) = -sum over faces containing edge(i,j) of cot(opposite angle)/2
%
% We use the negative sign convention so that L is positive semidefinite:
% off-diagonal entries are negative (attractive), diagonal entries positive.

% Edge (V2, V3): weight = cot1 (angle at V1)
% Edge (V1, V3): weight = cot2 (angle at V2)
% Edge (V1, V2): weight = cot3 (angle at V3)
% Each edge contributes to both L(i,j) and L(j,i) with the SAME weight.
ii = [Faces(:,2); Faces(:,3); Faces(:,1); Faces(:,3); Faces(:,1); Faces(:,2)];
jj = [Faces(:,3); Faces(:,2); Faces(:,3); Faces(:,1); Faces(:,2); Faces(:,1)];
ww = [cot1;       cot1;       cot2;       cot2;       cot3;       cot3      ];

% Off-diagonal: negative weights (L is positive semidefinite)
L = sparse(ii, jj, -ww, nV, nV);

% Diagonal: L(i,i) = -sum of off-diagonal row i
diagL = -full(sum(L, 2));
L = L + sparse(1:nV, 1:nV, diagL, nV, nV);

% Symmetrize (handles floating point asymmetry)
if Symmetrize
    L = (L + L') / 2;
end

%% ===== MASS MATRIX =====
% Triangle areas (half the cross product magnitude)
areas = 0.5 * dblArea;

switch MassType
    case 'barycentric'
        % Lumped mass: each vertex gets 1/3 of the area of each adjacent face.
        % M(i,i) = sum over faces containing i of area(face) / 3
        vertAreas = accumarray(Faces(:), repmat(areas, 3, 1) / 3, [nV, 1]);
        M = sparse(1:nV, 1:nV, vertAreas, nV, nV);

    case 'voronoi'
        % Voronoi dual cell areas. For non-obtuse triangles, the Voronoi
        % area at vertex i is:
        %   A_v(i) = (1/8) * sum_j [cot(alpha_ij) * ||e_ij||^2]
        % where the sum is over edges incident to i in each face.
        %
        % For obtuse triangles, we use Meyer et al. 2003 mixed Voronoi:
        %   - If angle at i is obtuse: A_v(i) = area/2
        %   - If another angle is obtuse: A_v(i) = area/4
        %   - Otherwise: standard Voronoi formula

        % Squared edge lengths
        l23sq = sum(e3.^2, 2);   % ||V3-V2||^2 = edge opposite V1
        l13sq = sum(e2.^2, 2);   % ||V3-V1||^2 = edge opposite V2
        l12sq = sum(e1.^2, 2);   % ||V2-V1||^2 = edge opposite V3

        % Angles at each vertex (using dot products; sign of dot = cos sign)
        % Obtuse if dot product < 0
        obtuse1 = (dot12 < 0);   % angle at V1 > 90 degrees
        obtuse2 = (dot_at2 < 0); % angle at V2 > 90 degrees
        obtuse3 = (dot_at3 < 0); % angle at V3 > 90 degrees
        anyObtuse = obtuse1 | obtuse2 | obtuse3;

        % Standard Voronoi areas (for non-obtuse faces):
        % Area at Vi = (1/8) * [cot(angle_j)*||e_ij||^2 + cot(angle_k)*||e_ik||^2]
        % But with our cot variables (already divided by dblArea):
        % Raw cotangents (without the 0.5 FEM scaling)
        raw_cot1 = dot12 ./ denom;
        raw_cot2 = dot_at2 ./ denom;
        raw_cot3 = dot_at3 ./ denom;

        % For vertex Vi, Voronoi area = (1/8) * sum over edges at Vi of
        % [cot(opposite angle) * ||edge||^2].
        %
        % V1: edges (V1,V2) and (V1,V3), opposite angles at V3 and V2
        %   Area_V1 = (1/8) * [cot(V3)*||V1-V2||^2 + cot(V2)*||V1-V3||^2]
        % V2: edges (V2,V1) and (V2,V3), opposite angles at V3 and V1
        %   Area_V2 = (1/8) * [cot(V3)*||V1-V2||^2 + cot(V1)*||V2-V3||^2]
        % V3: edges (V3,V1) and (V3,V2), opposite angles at V2 and V1
        %   Area_V3 = (1/8) * [cot(V2)*||V1-V3||^2 + cot(V1)*||V2-V3||^2]
        areaV1 = (1/8) * (raw_cot3 .* l12sq + raw_cot2 .* l13sq);
        areaV2 = (1/8) * (raw_cot3 .* l12sq + raw_cot1 .* l23sq);
        areaV3 = (1/8) * (raw_cot2 .* l13sq + raw_cot1 .* l23sq);

        % Mixed Voronoi: fix obtuse triangles (Meyer et al. 2003)
        % If the angle at the vertex is obtuse: area = total_area / 2
        % If another angle is obtuse: area = total_area / 4
        areaV1(obtuse1) = areas(obtuse1) / 2;
        areaV2(obtuse2) = areas(obtuse2) / 2;
        areaV3(obtuse3) = areas(obtuse3) / 2;

        % For vertices where ANOTHER angle is obtuse:
        otherObtuse1 = anyObtuse & ~obtuse1;
        otherObtuse2 = anyObtuse & ~obtuse2;
        otherObtuse3 = anyObtuse & ~obtuse3;
        areaV1(otherObtuse1) = areas(otherObtuse1) / 4;
        areaV2(otherObtuse2) = areas(otherObtuse2) / 4;
        areaV3(otherObtuse3) = areas(otherObtuse3) / 4;

        vertAreas = accumarray(Faces(:), [areaV1; areaV2; areaV3], [nV, 1]);
        M = sparse(1:nV, 1:nV, vertAreas, nV, nV);

    case 'galerkin'
        % Consistent FEM mass matrix (non-diagonal):
        %   M(i,i) += area(face) / 6   for each face containing vertex i
        %   M(i,j) += area(face) / 12  for each face containing edge (i,j)
        %
        % This gives M = integral of phi_i * phi_j over the mesh, where
        % phi are piecewise linear basis functions.
        ii_diag = Faces(:);
        jj_diag = Faces(:);
        ww_diag = repmat(areas / 6, 3, 1);

        ii_off = [Faces(:,1); Faces(:,2); Faces(:,3); Faces(:,2); Faces(:,3); Faces(:,1)];
        jj_off = [Faces(:,2); Faces(:,3); Faces(:,1); Faces(:,1); Faces(:,2); Faces(:,3)];
        ww_off = repmat(areas / 12, 6, 1);

        M = sparse([ii_diag; ii_off], [jj_diag; jj_off], [ww_diag; ww_off], nV, nV);

        if Symmetrize
            M = (M + M') / 2;
        end

    otherwise
        error('Unknown MassType: %s. Use ''barycentric'', ''voronoi'', or ''galerkin''.', MassType);
end

end
