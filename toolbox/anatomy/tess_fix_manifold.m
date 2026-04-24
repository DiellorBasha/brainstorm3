function [Vertices, Faces, nFixed, fixReport] = tess_fix_manifold(Vertices, Faces, varargin)
% TESS_FIX_MANIFOLD: Repair non-manifold defects in a triangulated surface.
%
% USAGE:  [Vertices, Faces, nFixed, fixReport] = tess_fix_manifold(Vertices, Faces)
%         [Vertices, Faces, nFixed, fixReport] = tess_fix_manifold(Vertices, Faces, 'Verbose', 1)
%
% DESCRIPTION:
%     Repairs non-manifold defects introduced by mesh downsampling (e.g.,
%     MATLAB's reducepatch). Handles:
%       1. Non-manifold edges (multiplicity >= 3): removes the spurious face
%          that causes the edge to be shared by too many triangles.
%       2. Resulting isolated/orphan faces and vertices are cleaned up.
%       3. Iterates until the mesh is edge-manifold or a maximum iteration
%          count is reached.
%
%     The repair strategy for non-manifold edges:
%       - For each edge shared by 3+ faces, identify which face to remove.
%       - Criterion: remove the face whose normal is most different from
%         the area-weighted average of its non-problematic neighbors.
%       - This preserves the local surface shape while eliminating the
%         topological defect.
%
%     Does NOT attempt to fill boundary holes (open edges). For eigenmode
%     computation on a cortical surface with two hemispheres, boundary edges
%     at the medial wall are expected and handled by the Laplacian assembly
%     with natural (Neumann) boundary conditions.
%
% INPUTS:
%     Vertices : [nVertices x 3] vertex positions
%     Faces    : [nFaces x 3] triangle vertex indices (1-based)
%
% OPTIONS:
%     Verbose  : (logical) Print progress. Default: 0
%     MaxIter  : (int) Maximum repair iterations. Default: 10
%
% OUTPUTS:
%     Vertices  : Repaired vertex array (vertices referenced by removed faces
%                 that become isolated are removed)
%     Faces     : Repaired face array
%     nFixed    : Number of faces removed
%     fixReport : Structure with repair details
%
% SEE ALSO: tess_check_manifold, tess_clean

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
Verbose = 0;
MaxIter = 10;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'verbose',  Verbose = varargin{i+1};
        case 'maxiter',  MaxIter = varargin{i+1};
    end
end

fixReport = struct();
fixReport.removedFaces = [];
fixReport.removedVertices = [];
fixReport.iterations = 0;
totalRemoved = 0;

%% ===== ITERATIVE REPAIR =====
for iter = 1:MaxIter
    nV = size(Vertices, 1);
    nF = size(Faces, 1);

    % Build edge incidence
    dE = [Faces(:,[1 2]); Faces(:,[2 3]); Faces(:,[3 1])];
    uE = sort(dE, 2);
    edgeKeys = uint64(uE(:,1)) * uint64(nV + 1) + uint64(uE(:,2));
    [~, ~, ic] = unique(edgeKeys);
    nE = max(ic);
    multiplicity = accumarray(ic, 1, [nE, 1]);

    % Find non-manifold edges
    nmEdgeIdx = find(multiplicity >= 3);
    if isempty(nmEdgeIdx)
        if Verbose
            fprintf('Iteration %d: No non-manifold edges. Done.\n', iter);
        end
        break;
    end

    if Verbose
        fprintf('Iteration %d: %d non-manifold edge(s), %d faces\n', iter, length(nmEdgeIdx), nF);
    end

    % Compute face normals and areas for scoring
    [faceNormals, faceAreas] = computeFaceNormalsAndAreas(Vertices, Faces);

    % Build face connectivity (which faces share edges)
    % For each directed edge occurrence, record its face index
    faceOfDE = mod((1:3*nF)' - 1, nF) + 1;

    % For each non-manifold edge, find adjacent faces and pick the one to remove
    facesToRemove = [];

    for k = 1:length(nmEdgeIdx)
        eIdx = nmEdgeIdx(k);
        % Find all directed edge occurrences belonging to this canonical edge
        occurrences = find(ic == eIdx);
        adjFaces = unique(faceOfDE(occurrences));

        if length(adjFaces) < 3
            continue;  % Already fixed by previous removal in this iteration
        end

        % Score each face: how well does its normal match the local neighborhood?
        % The "wrong" face typically has a normal that conflicts with its neighbors.
        worstScore = Inf;
        worstFace = adjFaces(1);

        for fi = 1:length(adjFaces)
            f = adjFaces(fi);
            fNorm = faceNormals(f, :);

            % Find neighbors of this face (faces sharing any edge, excluding
            % the other problematic faces at this non-manifold edge)
            fVerts = Faces(f, :);
            % Neighbor faces: share at least 2 vertices with f, but not in adjFaces
            neighbors = [];
            for vi = 1:3
                v = fVerts(vi);
                % Faces containing vertex v
                vFaces = find(any(Faces == v, 2))';
                neighbors = [neighbors, vFaces]; %#ok<AGROW>
            end
            neighbors = setdiff(unique(neighbors), adjFaces);

            if isempty(neighbors)
                % No clean neighbors — this face is likely the problem
                score = -Inf;
            else
                % Compute area-weighted average normal of clean neighbors
                neighborAreas = faceAreas(neighbors);
                avgNormal = sum(faceNormals(neighbors, :) .* neighborAreas, 1);
                nrm = norm(avgNormal);
                if nrm > 0
                    avgNormal = avgNormal / nrm;
                end
                % Score = dot product (higher = more consistent)
                score = dot(fNorm, avgNormal);
            end

            if score < worstScore
                worstScore = score;
                worstFace = f;
            end
        end

        facesToRemove(end+1) = worstFace; %#ok<AGROW>
    end

    % Remove the identified faces
    facesToRemove = unique(facesToRemove);
    if isempty(facesToRemove)
        break;
    end

    if Verbose
        fprintf('  Removing %d face(s).\n', length(facesToRemove));
    end

    fixReport.removedFaces = [fixReport.removedFaces; facesToRemove(:)];
    Faces(facesToRemove, :) = [];
    totalRemoved = totalRemoved + length(facesToRemove);

    % Remove isolated vertices (not referenced by any remaining face)
    usedVerts = unique(Faces(:));
    unusedVerts = setdiff(1:size(Vertices, 1), usedVerts);
    if ~isempty(unusedVerts)
        % Remap face indices
        vertMap = zeros(size(Vertices, 1), 1);
        vertMap(usedVerts) = 1:length(usedVerts);
        Vertices = Vertices(usedVerts, :);
        Faces = vertMap(Faces);
        fixReport.removedVertices = [fixReport.removedVertices; unusedVerts(:)];
        if Verbose
            fprintf('  Removed %d isolated vertex/vertices.\n', length(unusedVerts));
        end
    end

    fixReport.iterations = iter;
end

nFixed = totalRemoved;

%% ===== FIX ORIENTATION INCONSISTENCIES =====
% After removing non-manifold faces, some remaining faces may have
% inconsistent orientation (both faces traverse the shared edge in the
% same direction). Fix by flipping the winding of one face per bad edge.
%
% Algorithm: BFS-based orientation propagation.
%   1. Pick a seed face, mark its orientation as canonical.
%   2. For each neighbor sharing an interior edge, check if the shared
%      edge is traversed in opposite directions (consistent) or same
%      direction (inconsistent). If inconsistent, flip the neighbor.
%   3. Propagate through the entire connected component.

[Faces, nConsistencyFlips, nOutwardFlips] = fixOrientation(Vertices, Faces);
fixReport.nConsistencyFlips = nConsistencyFlips;
fixReport.nOutwardFlips = nOutwardFlips;
fixReport.nFlippedFaces = nConsistencyFlips + nOutwardFlips;
if Verbose
    if nConsistencyFlips > 0
        fprintf('Orientation: flipped %d face(s) for consistency.\n', nConsistencyFlips);
    end
    if nOutwardFlips > 0
        fprintf('Orientation: flipped %d face(s) for outward normals.\n', nOutwardFlips);
    end
end

%% ===== FINAL VALIDATION =====
[isManifold, finalReport] = tess_check_manifold(Vertices, Faces);
fixReport.finalValidation = finalReport;
fixReport.isManifold = isManifold;

if Verbose
    fprintf('\n--- Repair Summary ---\n');
    fprintf('  Removed %d face(s) in %d iteration(s).\n', nFixed, fixReport.iterations);
    fprintf('  Final: %d V, %d F\n', size(Vertices, 1), size(Faces, 1));
    if isManifold
        fprintf('  Result: MANIFOLD (edge-manifold, oriented, vertex-manifold)\n');
    else
        fprintf('  Result: Still has defects. Manual inspection needed.\n');
        for i = 1:length(finalReport.summary)
            fprintf('    %s\n', finalReport.summary{i});
        end
    end
    fprintf('----------------------\n');
end

end


%% ===== HELPER: FIX FACE ORIENTATION VIA BFS =====
function [F, nConsistencyFlips, nOutwardFlips] = fixOrientation(V, F)
    % BFS orientation propagation: make all faces globally consistent.
    %
    % After non-manifold edge removal, every undirected edge has at most
    % 2 adjacent faces. For a consistently-oriented pair, the directed edges
    % are (v1->v2) from face A and (v2->v1) from face B (twins). For an
    % inconsistently-oriented pair, both faces emit the same directed edge
    % (v1->v2).
    %
    % Strategy:
    %   1. Build face adjacency via undirected edges (canonical edge map).
    %   2. For each undirected edge shared by 2 faces, record whether the
    %      pair is consistent (twin) or inconsistent (same direction).
    %   3. BFS from a seed face. When traversing to a neighbor:
    %      - If the edge is "consistent" (twin): neighbor keeps its winding.
    %      - If the edge is "inconsistent" (same): flip the neighbor.
    %   4. After BFS, ensure outward orientation per component.

    nF = size(F, 1);
    nV = size(V, 1);

    % --- Build edge-face adjacency with orientation info ---
    % For each face, emit 3 directed edges: (v1,v2), (v2,v3), (v3,v1).
    v1s = [F(:,1); F(:,2); F(:,3)];
    v2s = [F(:,2); F(:,3); F(:,1)];
    faceIdx = repmat((1:nF)', 3, 1);  % which face this directed edge belongs to

    % Canonical (undirected) edge keys
    lo = min(v1s, v2s);
    hi = max(v1s, v2s);
    edgeKeys = uint64(lo) .* uint64(nV + 1) + uint64(hi);
    [~, ~, ic] = unique(edgeKeys);
    nE = max(ic);

    % For each canonical edge, collect adjacent faces and their directions.
    % Direction: +1 if the directed edge goes lo→hi, -1 if hi→lo.
    direction = ones(3*nF, 1);
    direction(v1s > v2s) = -1;

    % Build adjacency structures using accumarray.
    % For edges shared by exactly 2 faces, we need both face indices and
    % whether they have the same or opposite direction.
    %
    % Accumulate face indices per canonical edge (up to 2 per edge).
    multiplicity = accumarray(ic, 1, [nE, 1]);
    interiorEdges = find(multiplicity == 2);

    % For interior edges, find the two faces and their relative orientation.
    % Build sparse face-face adjacency with edge type (consistent/inconsistent).
    % adjType(f1, f2) = +1 if consistent (twin), -1 if inconsistent (same dir)
    adjI = zeros(0, 1);
    adjJ = zeros(0, 1);
    adjTypeVec = zeros(0, 1);

    if ~isempty(interiorEdges)
        % Sort directed edges by canonical edge index for grouping
        [icSorted, sortOrder] = sort(ic);
        faceIdxSorted = faceIdx(sortOrder);
        dirSorted = direction(sortOrder);

        % For interior edges (multiplicity==2), the two occurrences are
        % consecutive in the sorted list. Find the first occurrence of each.
        firstOcc = accumarray(icSorted, (1:length(icSorted))', [nE, 1], @min, 0);

        % Vectorized extraction: get both occurrences for all interior edges
        idx1 = firstOcc(interiorEdges);
        idx2 = idx1 + 1;

        f1 = faceIdxSorted(idx1);
        f2 = faceIdxSorted(idx2);
        d1 = dirSorted(idx1);
        d2 = dirSorted(idx2);

        % If d1 ~= d2 (opposite directions = twin edges): consistent (+1)
        % If d1 == d2 (same direction): inconsistent (-1)
        edgeTypes = ones(length(interiorEdges), 1);
        edgeTypes(d1 == d2) = -1;

        adjI = f1;
        adjJ = f2;
        adjTypeVec = edgeTypes;
    end

    % Build symmetric adjacency: adjType(f1,f2) = adjType(f2,f1)
    adjI_full = [adjI; adjJ];
    adjJ_full = [adjJ; adjI];
    adjType_full = [adjTypeVec; adjTypeVec];
    adjMatrix = sparse(adjI_full, adjJ_full, adjType_full, nF, nF);

    % --- BFS orientation propagation ---
    visited = false(nF, 1);
    flipped = false(nF, 1);
    queue = zeros(nF, 1);

    for seed = 1:nF
        if visited(seed)
            continue;
        end

        queue(1) = seed;
        qHead = 1;
        qTail = 1;
        visited(seed) = true;

        while qHead <= qTail
            curr = queue(qHead);
            qHead = qHead + 1;

            % Find all neighbors of curr in the adjacency
            [~, nbrs, types] = find(adjMatrix(curr, :));

            for ni = 1:length(nbrs)
                nbr = nbrs(ni);
                if visited(nbr)
                    continue;
                end

                edgeType = types(ni);

                % edgeType was computed BEFORE any flips. If curr was flipped,
                % the meaning of "consistent" is reversed for curr's edges.
                % If curr is flipped: consistent→inconsistent and vice versa.
                effectiveType = edgeType;
                if flipped(curr)
                    effectiveType = -effectiveType;
                end

                if effectiveType == -1
                    % Inconsistent with curr's current orientation → flip nbr
                    F(nbr, :) = F(nbr, [1 3 2]);
                    flipped(nbr) = true;
                end

                visited(nbr) = true;
                qTail = qTail + 1;
                queue(qTail) = nbr;
            end
        end
    end

    nConsistencyFlips = sum(flipped);
    nOutwardFlips = 0;

    % --- Ensure outward orientation per connected component ---
    % Use signed volume test: for an outward-oriented closed surface, the
    % signed volume V = (1/6) * sum_f (v1 . (v2 x v3)) > 0. For open
    % surfaces (e.g., cortical hemispheres with medial wall cut), the
    % signed volume is still a reliable majority indicator.
    % NOTE: The centroid-based dot-product test fails for highly convoluted
    % surfaces like cortex, where sulcal faces point toward the centroid.

    % Find connected components via vertex adjacency
    dE_new = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    uE_new = sort(dE_new, 2);
    Adj = sparse(double(uE_new(:,1)), double(uE_new(:,2)), 1, nV, nV);
    Adj = Adj + Adj';
    G = graph(Adj);
    compId = conncomp(G);

    % Map vertex components to face components
    faceComp = compId(F(:,1))';

    % Compute per-face signed volume contribution: (1/6) * v1 . (v2 x v3)
    v1 = V(F(:,1), :);
    v2 = V(F(:,2), :);
    v3 = V(F(:,3), :);
    signedVolPerFace = sum(v1 .* cross(v2, v3, 2), 2) / 6;

    for c = unique(faceComp)'
        mask = (faceComp == c);
        totalSignedVol = sum(signedVolPerFace(mask));
        if totalSignedVol < 0
            % Negative signed volume — normals point inward. Flip.
            idx = find(mask);
            F(idx, :) = F(idx, [1 3 2]);
            nOutwardFlips = nOutwardFlips + length(idx);
        end
    end
end


%% ===== HELPER: COMPUTE FACE NORMALS AND AREAS =====
function [normals, areas] = computeFaceNormalsAndAreas(V, F)
    v1 = V(F(:,2), :) - V(F(:,1), :);
    v2 = V(F(:,3), :) - V(F(:,1), :);
    cross_prod = [v1(:,2).*v2(:,3) - v1(:,3).*v2(:,2), ...
                  v1(:,3).*v2(:,1) - v1(:,1).*v2(:,3), ...
                  v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1)];
    areas = 0.5 * sqrt(sum(cross_prod.^2, 2));
    nrm = 2 * areas;  % = ||cross_prod||
    nrm(nrm == 0) = 1;  % avoid division by zero for degenerate faces
    normals = cross_prod ./ nrm;
end
