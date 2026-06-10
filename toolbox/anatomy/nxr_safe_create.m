function h = nxr_safe_create(V, F)
% NXR_SAFE_CREATE  Validate a mesh, then call nxr_compute('create', V, F).
%
% USAGE:  h = nxr_safe_create(V, F)
%
% geometry-central's ManifoldSurfaceMesh constructor SEGFAULTS (taking down
% the whole MATLAB process) when handed malformed mesh data: out-of-range or
% non-positive face indices, degenerate faces, non-manifold edges, or
% inconsistent winding. This wrapper validates V/F in MATLAB FIRST so a bad
% mesh raises a clean, catchable error instead of crashing the session.
%
% Use this everywhere nxr_compute('create', ...) would be called on a mesh
% that is not already known-good.
%
% INPUT:
%   V : [nV x 3] real, finite vertex coordinates
%   F : [nF x 3] triangle indices, 1-based, in [1, nV]
%
% OUTPUT:
%   h : opaque nxr handle (caller must nxr_compute('destroy', h))

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% =============================================================================@

    % ---- vertices ----
    if ~ismatrix(V) || size(V,2) ~= 3 || isempty(V)
        error('nxr_safe_create:vertices', 'V must be a non-empty [nV x 3] matrix.');
    end
    if ~isreal(V) || ~all(isfinite(V(:)))
        error('nxr_safe_create:vertices', 'V must be real and finite (no NaN/Inf).');
    end
    nV = size(V, 1);

    % ---- faces ----
    if ~ismatrix(F) || size(F,2) ~= 3 || isempty(F)
        error('nxr_safe_create:faces', 'F must be a non-empty [nF x 3] matrix.');
    end
    if ~all(isfinite(F(:))) || any(F(:) ~= floor(F(:)))
        error('nxr_safe_create:faces', 'F must contain only integer indices.');
    end
    if any(F(:) < 1)
        error('nxr_safe_create:faceBase', ...
            ['F has indices < 1 (got min=%d). nxr expects 1-based face indices; ', ...
             'a 0-based mesh underflows inside geometry-central and segfaults.'], min(F(:)));
    end
    if any(F(:) > nV)
        error('nxr_safe_create:faceRange', ...
            'F references vertex %d but V has only %d vertices.', max(F(:)), nV);
    end

    % ---- degenerate faces (repeated vertex within a triangle) ----
    degen = (F(:,1) == F(:,2)) | (F(:,2) == F(:,3)) | (F(:,1) == F(:,3));
    if any(degen)
        error('nxr_safe_create:degenerate', ...
            '%d face(s) are degenerate (a vertex repeated within the triangle).', nnz(degen));
    end

    % ---- edge-manifold check (every undirected edge shared by <= 2 faces) ----
    dirEdges = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];      % directed half-edges
    undEdges = sort(dirEdges, 2);                          % undirected edges
    [~, ~, iu] = unique(undEdges, 'rows');
    if any(accumarray(iu, 1) > 2)
        error('nxr_safe_create:nonManifold', ...
            ['Mesh is non-manifold: an edge is shared by more than two faces. ', ...
             'geometry-central''s ManifoldSurfaceMesh cannot represent this.']);
    end

    % ---- orientation check (no directed half-edge appears twice) ----
    [~, ~, id] = unique(dirEdges, 'rows');
    if any(accumarray(id, 1) > 1)
        error('nxr_safe_create:orientation', ...
            ['Mesh has inconsistent face orientation (a directed edge appears ', ...
             'in two faces with the same winding). Faces must be consistently oriented.']);
    end

    % All checks passed — safe to construct the manifold.
    h = nxr_compute('create', V, F);
end
