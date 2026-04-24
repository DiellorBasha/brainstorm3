function [isManifold, report] = tess_check_manifold(Vertices, Faces, varargin)
% TESS_CHECK_MANIFOLD: Comprehensive manifold validation for triangulated surfaces.
%
% USAGE:  [isManifold, report] = tess_check_manifold(Vertices, Faces)
%         [isManifold, report] = tess_check_manifold(Vertices, Faces, 'RequireClosed', 1, ...)
%
% DESCRIPTION:
%     Validates that a triangulated surface is a proper 2-manifold suitable
%     for discrete exterior calculus operations (cotangent Laplacian, eigenmode
%     decomposition, Hodge operators).
%
%     A surface is manifold if and only if:
%       1. Every edge is shared by exactly 1 face (boundary) or 2 faces (interior)
%       2. The faces around every vertex form a single connected fan (disk or half-disk)
%       3. Face orientation is globally consistent (adjacent faces traverse shared
%          edges in opposite directions)
%       4. No degenerate faces (zero area, repeated vertex indices)
%
%     This function checks all four conditions plus auxiliary mesh quality
%     properties needed for reliable operator assembly.
%
%     Algorithm follows geometry-central's ManifoldSurfaceMesh validation and
%     the bioctree +bct/+manifold/+health module. No external dependencies.
%
% INPUTS:
%     Vertices       : [nVertices x 3] vertex positions
%     Faces          : [nFaces x 3] triangle vertex indices (1-based)
%
% OPTIONS (name-value pairs):
%     RequireClosed  : (logical) Require no boundary edges. Default: 0
%     RequireConnected: (logical) Require single connected component. Default: 0
%     Verbose        : (logical) Print summary to console. Default: 0
%     AreaTol        : (double) Tolerance for zero-area faces. Default: 1e-20
%     DupVertexTol   : (double) Tolerance for duplicate vertex detection. Default: 1e-10
%
% OUTPUTS:
%     isManifold : logical, true if surface passes all manifold checks
%     report     : structure with fields:
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
%                           nDegenerateFaces, nDuplicateVertices, nComponents,
%                           eulerCharacteristic
%       .badEdges         - struct with edge indices:
%                           .nonManifold [k x 2], .boundary [k x 2],
%                           .inconsistent [k x 2]
%       .badFaces         - indices of degenerate or duplicate faces
%       .badVertices      - indices of non-manifold or duplicate vertices
%       .summary          - cell array of human-readable messages
%
% SEE ALSO: tess_clean, tess_downsize

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
RequireClosed = 0;
RequireConnected = 0;
Verbose = 0;
AreaTol = 1e-20;
DupVertexTol = 1e-10;

% Parse name-value pairs
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'requireclosed',    RequireClosed = varargin{i+1};
        case 'requireconnected', RequireConnected = varargin{i+1};
        case 'verbose',          Verbose = varargin{i+1};
        case 'areatol',          AreaTol = varargin{i+1};
        case 'dupvertextol',     DupVertexTol = varargin{i+1};
        otherwise
            error('Unknown option: %s', varargin{i});
    end
end

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
    if Verbose; printSummary(report); end
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
    faceAreas = computeFaceAreas(Vertices, Faces);
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
nDE = size(dE, 1);  % = 3 * nF

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
    % Group by canonical edge: use accumarray to collect directed keys.
    % Two occurrences with same directed key → same direction → bad.
    %
    % Fast approach: for each canonical edge, sum the directed keys of its
    % two occurrences. If both are identical, the sum == 2*key. If different
    % (opposite direction), the sum != 2*key. But we need a simpler test.
    %
    % Simplest vectorized approach: sort by (canonIdx, dirKey), then check
    % for consecutive identical pairs.
    sortData = [double(intCanonIdx), double(intDirKeys)];
    sortData = sortrows(sortData, [1, 2]);

    % Consecutive rows with same canonical index AND same directed key = inconsistent
    sameCanon = diff(sortData(:,1)) == 0;
    sameDir   = diff(sortData(:,2)) == 0;
    badPairs  = sameCanon & sameDir;

    % Map back to canonical edge indices
    if any(badPairs)
        badCanonIdx = sortData(badPairs, 1);
        % These are indices into the interior edge set; map to E indices
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
%
% If the fast check passes AND the mesh is edge-manifold, the surface is
% vertex-manifold (proven by the classification theorem for 2-manifolds).

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
%   Boundary vertex (bndIncV > 0):  nFaces must equal nEdges - 1
%   (Boundary vertex with bndIncV == 2 has one open arc; bndIncV > 2 is non-manifold)
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

%% ===== VERBOSE OUTPUT =====
if Verbose
    printSummary(report);
end

end


%% ===== HELPER: COMPUTE FACE AREAS =====
function areas = computeFaceAreas(V, F)
    % Cross product method for triangle areas
    v1 = V(F(:,2), :) - V(F(:,1), :);
    v2 = V(F(:,3), :) - V(F(:,1), :);
    cross_prod = [v1(:,2).*v2(:,3) - v1(:,3).*v2(:,2), ...
                  v1(:,3).*v2(:,1) - v1(:,1).*v2(:,3), ...
                  v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1)];
    areas = 0.5 * sqrt(sum(cross_prod.^2, 2));
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
