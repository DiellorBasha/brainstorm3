function [Vertices, Faces, isManifold, report] = tess_repair(Vertices, Faces, varargin)
% TESS_REPAIR: Validate and (optionally) repair a triangulated surface as a 2-manifold.
%
% USAGE:  [~, ~, isManifold, report]            = tess_repair(Vertices, Faces)              % check only
%         [Vertices, Faces, isManifold, report] = tess_repair(Vertices, Faces, 'Repair', 1) % check + repair
%         [Vertices, Faces, isManifold, report] = tess_repair(Vertices, Faces, 'RequireClosed', 1, ...)
%
% DESCRIPTION:
%     Single entry point for 2-manifold validation and repair of a triangulated
%     surface (e.g. a cortical surface produced by downsampling). Replaces the
%     former tess_check_manifold (validation) and tess_fix_manifold (repair).
%
%     By DEFAULT the function only VALIDATES: it checks that the surface is a
%     proper 2-manifold suitable for discrete exterior calculus operations
%     (cotangent Laplacian, eigenmode decomposition, Hodge operators) and
%     returns the verdict without touching the mesh.
%
%     A surface is manifold if and only if:
%       1. Every edge is shared by exactly 1 face (boundary) or 2 faces (interior)
%       2. The faces around every vertex form a single connected fan (disk or half-disk)
%       3. Face orientation is globally consistent (adjacent faces traverse shared
%          edges in opposite directions)
%       4. No degenerate faces (zero area, repeated vertex indices)
%
%     The validation checks all four conditions plus auxiliary mesh quality
%     properties needed for reliable operator assembly. The algorithm follows
%     geometry-central's ManifoldSurfaceMesh validation and the bioctree
%     +bct/+manifold/+health module. No external dependencies.
%
%     With 'Repair', 1 the function validates first and, ONLY IF the surface is
%     defective, repairs the non-manifold defects introduced by mesh
%     downsampling (e.g. MATLAB's reducepatch), then re-validates and returns
%     the fixed mesh. An already-manifold surface is returned unchanged (no-op).
%     Repair handles:
%       1. Non-manifold edges (multiplicity >= 3): removes the spurious face
%          that causes the edge to be shared by too many triangles. Criterion:
%          remove the face whose normal is most different from the area-weighted
%          average of its non-problematic neighbors (preserves local shape).
%       2. Resulting isolated/orphan faces and vertices are cleaned up.
%       3. Iterates until edge-manifold or MaxIter is reached.
%       4. Orientation: BFS propagation makes all faces consistently oriented,
%          then a signed-volume test ensures outward normals per component.
%     Repair does NOT fill boundary holes (open edges). For eigenmode
%     computation on a cortical surface with two hemispheres, boundary edges at
%     the medial wall are expected and handled by the Laplacian assembly with
%     natural (Neumann) boundary conditions.
%
% INPUTS:
%     Vertices : [nVertices x 3] vertex positions
%     Faces    : [nFaces x 3] triangle vertex indices (1-based)
%
% OPTIONS (name-value pairs):
%     Repair          : (logical) Repair defects if found. Default: 0 (check only)
%     RequireClosed   : (logical) Treat boundary edges as a failure. Default: 0
%     RequireConnected: (logical) Require single connected component. Default: 0
%     Verbose         : (logical) Print summary to console. Default: 0
%     AreaTol         : (double) Tolerance for zero-area faces. Default: 1e-20
%     DupVertexTol    : (double) Tolerance for duplicate vertex detection. Default: 1e-10
%     MaxIter         : (int) Maximum repair iterations. Default: 10
%
% OUTPUTS:
%     Vertices   : Vertex array. Unchanged unless 'Repair' removed isolated vertices.
%     Faces      : Face array. Unchanged unless 'Repair' removed/reoriented faces.
%     isManifold : logical, true if the RETURNED surface passes all manifold checks
%                  (post-repair state when 'Repair' was requested).
%     report     : structure describing the RETURNED mesh, with fields:
%       .ok               - overall pass/fail (same as isManifold)
%       .isEdgeManifold   - no edge shared by 3+ faces
%       .isOriented       - consistent face orientation
%       .isVertexManifold - disk topology at every vertex
%       .isClosed         - no boundary edges (informational unless RequireClosed)
%       .isConnected      - single connected component
%       .noDegenerateFaces- no zero-area or repeated-vertex faces
%       .noDuplicateFaces - no identical face triangles
%       .noDuplicateVerts - no coincident vertices
%       .stats            - struct with mesh statistics:
%                           nVertices, nFaces, nEdges, nBoundaryEdges,
%                           nNonManifoldEdges, nInconsistentEdges,
%                           nDegenerateFaces, nDuplicateFaces, nDuplicateVertices,
%                           nComponents, eulerCharacteristic
%       .badEdges         - struct with edge vertex pairs:
%                           .nonManifold [k x 2], .boundary [k x 2],
%                           .inconsistent [k x 2]
%       .badFaces         - indices of degenerate or duplicate faces
%       .badVertices      - indices of non-manifold or duplicate vertices
%       .summary          - cell array of human-readable messages
%       .repair           - (only when 'Repair' was requested) struct with:
%                           .performed (logical, false if mesh was already clean),
%                           .nFixed, .iterations, .removedFaces, .removedVertices,
%                           .nConsistencyFlips, .nOutwardFlips, .nFlippedFaces,
%                           .validationBefore (the pre-repair validation report)
%
% SEE ALSO: tess_laplacian, tess_eigenmodes, tess_clean, tess_downsize

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
% Defaults
Repair = 0;
RequireClosed = 0;
RequireConnected = 0;
Verbose = 0;
AreaTol = 1e-20;
DupVertexTol = 1e-10;
MaxIter = 10;

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'repair',           Repair = varargin{i+1};
        case 'requireclosed',    RequireClosed = varargin{i+1};
        case 'requireconnected', RequireConnected = varargin{i+1};
        case 'verbose',          Verbose = varargin{i+1};
        case 'areatol',          AreaTol = varargin{i+1};
        case 'dupvertextol',     DupVertexTol = varargin{i+1};
        case 'maxiter',          MaxIter = varargin{i+1};
        otherwise
            error('Unknown option: %s', varargin{i});
    end
end

%% ===== VALIDATE CURRENT MESH =====
[isManifold, report] = doCheck(Vertices, Faces, RequireClosed, RequireConnected, AreaTol, DupVertexTol);

%% ===== REPAIR IF REQUESTED =====
if Repair && ~isManifold
    validationBefore = report;
    if Verbose
        fprintf('BST> tess_repair: Mesh has defects, repairing...\n');
    end
    [Vertices, Faces, repairDetails] = doRepair(Vertices, Faces, MaxIter, Verbose);
    % Re-validate the repaired mesh: the report must describe what the caller gets back.
    [isManifold, report] = doCheck(Vertices, Faces, RequireClosed, RequireConnected, AreaTol, DupVertexTol);
    repairDetails.performed = true;
    repairDetails.validationBefore = validationBefore;
    report.repair = repairDetails;
elseif Repair
    % Already manifold: repair is a no-op, mesh returned unchanged.
    report.repair = struct(...
        'performed',         false, ...
        'nFixed',            0, ...
        'iterations',        0, ...
        'removedFaces',      [], ...
        'removedVertices',   [], ...
        'nConsistencyFlips', 0, ...
        'nOutwardFlips',     0, ...
        'nFlippedFaces',     0, ...
        'validationBefore',  report);
end

%% ===== VERBOSE OUTPUT =====
if Verbose
    printSummary(report);
    if isfield(report, 'repair') && report.repair.performed
        fprintf('BST> tess_repair: Repaired %d face(s) in %d iteration(s), %d orientation flip(s). Final: %d V, %d F.\n', ...
            report.repair.nFixed, report.repair.iterations, report.repair.nFlippedFaces, ...
            size(Vertices, 1), size(Faces, 1));
    end
end

end


%% ========================================================================
%  VALIDATION
%  ========================================================================
function [isManifold, report] = doCheck(Vertices, Faces, RequireClosed, RequireConnected, AreaTol, DupVertexTol)
% Comprehensive 2-manifold validation (no console output, no mesh modification).

%% ===== INPUT VALIDATION =====
if size(Vertices, 2) ~= 3
    error('Vertices must be [N x 3].');
end
if size(Faces, 2) ~= 3
    error('Faces must be [N x 3].');
end

nV = size(Vertices, 1);
nF = size(Faces, 1);

% Initialize report
report = struct();
report.summary = {};

%% ===== CHECK 1: FACE INDEX VALIDITY =====
% All face indices must be in [1, nV] and integer-valued
facesValid = all(Faces(:) >= 1) && all(Faces(:) <= nV) && ...
             all(Faces(:) == floor(Faces(:)));
if ~facesValid
    report.ok = false;
    report.summary{end+1} = sprintf('ERROR: Face indices out of range [1, %d] or non-integer.', nV);
    report = fillDefaults(report, nV, nF);
    isManifold = false;
    return;
end

%% ===== CHECK 2: DEGENERATE FACES =====
% A face is degenerate if any two vertex indices are identical
degMask = (Faces(:,1) == Faces(:,2)) | ...
          (Faces(:,2) == Faces(:,3)) | ...
          (Faces(:,3) == Faces(:,1));
degFaceIdx = find(degMask);
nDegFaces = length(degFaceIdx);
report.noDegenerateFaces = (nDegFaces == 0);
if nDegFaces > 0
    report.summary{end+1} = sprintf('ERROR: %d degenerate face(s) (repeated vertex indices).', nDegFaces);
end

% Also check for near-zero-area faces (requires vertices)
if nDegFaces == 0
    [~, faceAreas] = computeFaceNormalsAndAreas(Vertices, Faces);
    zeroAreaIdx = find(faceAreas < AreaTol);
    if ~isempty(zeroAreaIdx)
        degFaceIdx = zeroAreaIdx;
        nDegFaces = length(zeroAreaIdx);
        report.noDegenerateFaces = false;
        report.summary{end+1} = sprintf('WARNING: %d near-zero-area face(s) (area < %.1e).', nDegFaces, AreaTol);
    end
end

%% ===== CHECK 3: DUPLICATE FACES =====
sortedF = sort(Faces, 2);
[~, uniqueIdx] = unique(sortedF, 'rows', 'first');
dupFaceIdx = setdiff(1:nF, uniqueIdx);
nDupFaces = length(dupFaceIdx);
report.noDuplicateFaces = (nDupFaces == 0);
if nDupFaces > 0
    report.summary{end+1} = sprintf('ERROR: %d duplicate face(s).', nDupFaces);
end

%% ===== BUILD EDGE INCIDENCE =====
% Directed edges from faces: edges (1→2), (2→3), (3→1) per face
dE = [Faces(:,[1 2]); Faces(:,[2 3]); Faces(:,[3 1])];  % [3*nF x 2]

% Canonical (undirected) edge: always store (min, max)
uE = sort(dE, 2);  % [3*nF x 2]

% Unique edges and their mapping
% Encode edge as single integer for fast grouping: key = v1 * (nV+1) + v2
edgeKeys = uint64(uE(:,1)) * uint64(nV + 1) + uint64(uE(:,2));
[uniqueKeys, ~, ic] = unique(edgeKeys);
nE = length(uniqueKeys);

% Recover canonical edge vertex pairs
canonV1 = floor(double(uniqueKeys) / double(nV + 1));
canonV2 = mod(double(uniqueKeys), double(nV + 1));
E = [canonV1(:), canonV2(:)];

% Multiplicity: how many faces share each edge
multiplicity = accumarray(ic, 1, [nE, 1]);

%% ===== CHECK 4: EDGE MANIFOLDNESS =====
% Interior edge: multiplicity == 2, Boundary edge: multiplicity == 1
% Non-manifold: multiplicity >= 3
nonManifoldIdx = find(multiplicity >= 3);
nNonManifold = length(nonManifoldIdx);
report.isEdgeManifold = (nNonManifold == 0);
if nNonManifold > 0
    report.summary{end+1} = sprintf('ERROR: %d non-manifold edge(s) (shared by 3+ faces).', nNonManifold);
end

%% ===== CHECK 5: BOUNDARY EDGES =====
boundaryIdx = find(multiplicity == 1);
nBoundary = length(boundaryIdx);
report.isClosed = (nBoundary == 0);
if nBoundary > 0
    msg = sprintf('INFO: %d boundary edge(s) (open mesh).', nBoundary);
    if RequireClosed
        msg = strrep(msg, 'INFO', 'ERROR');
    end
    report.summary{end+1} = msg;
end

%% ===== CHECK 6: ORIENTATION CONSISTENCY =====
% For each interior edge (multiplicity == 2), the two directed edges from
% the two adjacent faces should traverse in opposite directions:
%   Face A contributes [u, v], Face B should contribute [v, u].
% If both contribute [u, v], orientation is inconsistent.
%
% Vectorized algorithm: encode directed edges as keys, group by canonical
% edge, and check for duplicate directed keys within each group.

interiorEdgeIdx = find(multiplicity == 2);
nInconsistent = 0;
inconsistentEdgeIdx = [];

if ~isempty(interiorEdgeIdx)
    % Encode directed edges as unique keys (order-dependent, unlike canonical)
    dirKeys = uint64(dE(:,1)) * uint64(nV + 1) + uint64(dE(:,2));

    % For each directed edge occurrence, get its canonical edge index (ic)
    % and whether it's an interior edge
    isInterior = (multiplicity(ic) == 2);

    % Keep only directed edges that belong to interior canonical edges
    intDirKeys = dirKeys(isInterior);
    intCanonIdx = ic(isInterior);

    % For interior edges with multiplicity 2: each canonical edge has exactly
    % 2 directed edge occurrences. If both have the same directed key, they
    % traverse in the same direction → inconsistent.
    %
    % Sort by (canonIdx, dirKey), then check for consecutive identical pairs.
    sortData = [double(intCanonIdx), double(intDirKeys)];
    sortData = sortrows(sortData, [1, 2]);

    % Consecutive rows with same canonical index AND same directed key = inconsistent
    sameCanon = diff(sortData(:,1)) == 0;
    sameDir   = diff(sortData(:,2)) == 0;
    badPairs  = sameCanon & sameDir;

    % Map back to canonical edge indices
    if any(badPairs)
        badCanonIdx = sortData(badPairs, 1);
        inconsistentEdgeIdx = unique(badCanonIdx);
        nInconsistent = length(inconsistentEdgeIdx);
    end
end

report.isOriented = (nInconsistent == 0);
if nInconsistent > 0
    report.summary{end+1} = sprintf('ERROR: %d edge(s) with inconsistent face orientation.', nInconsistent);
end

%% ===== CHECK 7: VERTEX MANIFOLDNESS =====
% Each vertex must have a single connected fan of faces around it.
%
% Fast necessary condition (O(V) with precomputed data):
%   For an interior vertex: nIncidentFaces == nIncidentEdges (closed fan)
%   For a boundary vertex:  nIncidentFaces == nIncidentEdges - 1 (open arc)
%   Any other count means disconnected fan → non-manifold vertex.
%
% This catches "bowtie" / "pinch point" vertices where two surface sheets
% meet at a single vertex. It requires edge manifoldness to already be checked.

report.isVertexManifold = true;
badVertexIdx = [];

% Count incident edges per vertex from canonical edge list
edgeIncV = accumarray([E(:,1); E(:,2)], 1, [nV, 1]);

% Count incident faces per vertex
faceIncV = accumarray([Faces(:,1); Faces(:,2); Faces(:,3)], 1, [nV, 1]);

% Count boundary edges per vertex (edges with multiplicity == 1)
if nBoundary > 0
    bndEdges = E(boundaryIdx, :);
    bndIncV = accumarray([bndEdges(:,1); bndEdges(:,2)], 1, [nV, 1]);
else
    bndIncV = zeros(nV, 1);
end

% For each vertex, check the fan condition:
%   Interior vertex (bndIncV == 0): nFaces must equal nEdges
%   Boundary vertex (bndIncV == 2): nFaces must equal nEdges - 1
%   (bndIncV == 1 or bndIncV >= 3 is non-manifold)
for v = 1:nV
    nf = faceIncV(v);
    ne = edgeIncV(v);
    nb = bndIncV(v);

    if nf == 0 && ne == 0
        continue;  % Isolated vertex — not a topology error
    end

    if nb == 0
        % Interior vertex: closed fan requires nFaces == nEdges
        if nf ~= ne
            badVertexIdx(end+1) = v; %#ok<AGROW>
        end
    elseif nb == 2
        % Standard boundary vertex (single arc): nFaces == nEdges - 1
        if nf ~= ne - 1
            badVertexIdx(end+1) = v; %#ok<AGROW>
        end
    else
        % nb == 1 or nb >= 3: topologically inconsistent for a manifold vertex
        badVertexIdx(end+1) = v; %#ok<AGROW>
    end
end

if ~isempty(badVertexIdx)
    report.isVertexManifold = false;
    report.summary{end+1} = sprintf('ERROR: %d non-manifold vertex/vertices (disconnected fan).', length(badVertexIdx));
end

%% ===== CHECK 8: DUPLICATE VERTICES =====
nDupVerts = 0;
dupVertIdx = [];
if DupVertexTol > 0
    % Round to tolerance and find duplicates
    Vround = round(Vertices / DupVertexTol) * DupVertexTol;
    [~, uniqueVIdx] = unique(Vround, 'rows', 'first');
    dupVertIdx = setdiff(1:nV, uniqueVIdx);
    nDupVerts = length(dupVertIdx);
end
report.noDuplicateVerts = (nDupVerts == 0);
if nDupVerts > 0
    report.summary{end+1} = sprintf('WARNING: %d duplicate vertex/vertices (within tol %.1e).', nDupVerts, DupVertexTol);
end

%% ===== CHECK 9: CONNECTIVITY =====
% Build symmetric adjacency matrix from edges and use MATLAB's graph/conncomp
A = sparse(double(E(:,1)), double(E(:,2)), 1, nV, nV);
A = A + A';  % symmetrize
G = graph(A);
componentId = conncomp(G);
nComponents = max(componentId);

report.isConnected = (nComponents == 1);
if nComponents > 1
    msg = sprintf('INFO: %d connected component(s).', nComponents);
    if RequireConnected
        msg = strrep(msg, 'INFO', 'ERROR');
    end
    report.summary{end+1} = msg;
end

%% ===== CHECK 10: EULER CHARACTERISTIC =====
% For a closed orientable manifold: V - E + F = 2 * (1 - genus)
% For a sphere (genus 0): V - E + F = 2
eulerChar = nV - nE + nF;
report.summary{end+1} = sprintf('INFO: Euler characteristic = %d (V=%d, E=%d, F=%d).', eulerChar, nV, nE, nF);

%% ===== AGGREGATE RESULTS =====
% A surface is manifold if edges are manifold, orientation is consistent,
% vertices are manifold, and there are no degenerate/duplicate faces.
isManifold = report.isEdgeManifold && ...
             report.isOriented && ...
             report.isVertexManifold && ...
             report.noDegenerateFaces && ...
             report.noDuplicateFaces;

% Also fail if closed/connected required and not met
if RequireClosed && ~report.isClosed
    isManifold = false;
end
if RequireConnected && ~report.isConnected
    isManifold = false;
end

report.ok = isManifold;

% Statistics
report.stats = struct(...
    'nVertices',          nV, ...
    'nFaces',             nF, ...
    'nEdges',             nE, ...
    'nBoundaryEdges',     nBoundary, ...
    'nNonManifoldEdges',  nNonManifold, ...
    'nInconsistentEdges', nInconsistent, ...
    'nDegenerateFaces',   nDegFaces, ...
    'nDuplicateFaces',    nDupFaces, ...
    'nDuplicateVertices', nDupVerts, ...
    'nComponents',        nComponents, ...
    'eulerCharacteristic', eulerChar);

% Bad element indices
report.badEdges = struct(...
    'nonManifold',   E(nonManifoldIdx, :), ...
    'boundary',      E(boundaryIdx, :), ...
    'inconsistent',  E(inconsistentEdgeIdx, :));
report.badFaces = union(degFaceIdx(:)', dupFaceIdx(:)');
report.badVertices = union(badVertexIdx(:)', dupVertIdx(:)');

% Final summary line
if isManifold
    report.summary{end+1} = 'PASS: Surface is a valid 2-manifold.';
else
    report.summary{end+1} = 'FAIL: Surface is NOT a valid 2-manifold.';
end

end


%% ========================================================================
%  REPAIR
%  ========================================================================
function [Vertices, Faces, repairDetails] = doRepair(Vertices, Faces, MaxIter, Verbose)
% Repair non-manifold defects: remove spurious faces, clean orphans, fix
% orientation. Caller re-validates afterwards.

repairDetails = struct();
repairDetails.removedFaces = [];
repairDetails.removedVertices = [];
repairDetails.iterations = 0;
totalRemoved = 0;

%% ===== ITERATIVE EDGE REPAIR =====
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
            fprintf('  Iteration %d: No non-manifold edges. Done.\n', iter);
        end
        break;
    end

    if Verbose
        fprintf('  Iteration %d: %d non-manifold edge(s), %d faces\n', iter, length(nmEdgeIdx), nF);
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
        fprintf('    Removing %d face(s).\n', length(facesToRemove));
    end

    repairDetails.removedFaces = [repairDetails.removedFaces; facesToRemove(:)];
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
        repairDetails.removedVertices = [repairDetails.removedVertices; unusedVerts(:)];
        if Verbose
            fprintf('    Removed %d isolated vertex/vertices.\n', length(unusedVerts));
        end
    end

    repairDetails.iterations = iter;
end

repairDetails.nFixed = totalRemoved;

%% ===== FIX ORIENTATION INCONSISTENCIES =====
% After removing non-manifold faces, some remaining faces may have
% inconsistent orientation (both faces traverse the shared edge in the
% same direction). Fix by flipping the winding of one face per bad edge,
% then ensure outward normals per connected component.

[Faces, nConsistencyFlips, nOutwardFlips] = fixOrientation(Vertices, Faces);
repairDetails.nConsistencyFlips = nConsistencyFlips;
repairDetails.nOutwardFlips = nOutwardFlips;
repairDetails.nFlippedFaces = nConsistencyFlips + nOutwardFlips;
if Verbose
    if nConsistencyFlips > 0
        fprintf('  Orientation: flipped %d face(s) for consistency.\n', nConsistencyFlips);
    end
    if nOutwardFlips > 0
        fprintf('  Orientation: flipped %d face(s) for outward normals.\n', nOutwardFlips);
    end
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
    multiplicity = accumarray(ic, 1, [nE, 1]);
    interiorEdges = find(multiplicity == 2);

    % For interior edges, find the two faces and their relative orientation.
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


%% ===== HELPER: FILL DEFAULTS FOR EARLY EXIT =====
function report = fillDefaults(report, nV, nF)
    report.ok = false;
    report.isEdgeManifold = false;
    report.isOriented = false;
    report.isVertexManifold = false;
    report.isClosed = false;
    report.isConnected = false;
    report.noDegenerateFaces = false;
    report.noDuplicateFaces = false;
    report.noDuplicateVerts = false;
    report.stats = struct('nVertices', nV, 'nFaces', nF, 'nEdges', 0, ...
        'nBoundaryEdges', 0, 'nNonManifoldEdges', 0, 'nInconsistentEdges', 0, ...
        'nDegenerateFaces', 0, 'nDuplicateFaces', 0, 'nDuplicateVertices', 0, ...
        'nComponents', 0, 'eulerCharacteristic', 0);
    report.badEdges = struct('nonManifold', [], 'boundary', [], 'inconsistent', []);
    report.badFaces = [];
    report.badVertices = [];
end


%% ===== HELPER: PRINT SUMMARY =====
function printSummary(report)
    fprintf('\n--- Manifold Validation Report ---\n');
    for i = 1:length(report.summary)
        fprintf('  %s\n', report.summary{i});
    end
    fprintf('----------------------------------\n\n');
end
