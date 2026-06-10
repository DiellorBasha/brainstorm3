function [SurfaceFile, Vertices, Faces] = bst_canonical_cortex(nVertTarget)
% BST_CANONICAL_CORTEX  The canonical cortical-surface test manifold.
%
% USAGE:  [SurfaceFile, Vertices, Faces] = bst_canonical_cortex(nVertTarget)
%         SurfaceFile = bst_canonical_cortex()        % defaults to 20484
%
% This is the SINGLE source of truth for the cortical manifold used in tests
% and probes. It loads the real Cortex surface from the currently-loaded
% Brainstorm protocol (starting Brainstorm headless if needed).
%
% HARD RULE — never hand-build a mesh for nxr_compute('create'):
%   Any test or probe that calls nxr_compute('create', V, F) MUST obtain its
%   mesh from here (then split per-hemisphere with tess_hemisplit and reindex
%   locally, exactly as tess_manifold does). Hand-built meshes (octahedron,
%   sphere, tetrahedron, ...) are BANNED as input to nxr_compute('create'):
%   geometry-central's ManifoldSurfaceMesh constructor SEGFAULTS (crashing
%   MATLAB) on malformed input rather than throwing. Synthetic data is allowed
%   ONLY for pure-algorithm tests that never call nxr_compute('create').
%
% OUTPUT:
%   SurfaceFile : protocol-relative path to the canonical Cortex surface
%   Vertices    : [nVert x 3] vertex coordinates (loaded from the surface)
%   Faces       : [nFace x 3] 1-based triangle indices (loaded from the surface)

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% =============================================================================@

    if nargin < 1 || isempty(nVertTarget)
        nVertTarget = 20484;
    end

    % Ensure a protocol is loaded
    if ~brainstorm('status')
        brainstorm nogui;
    end

    SurfaceFile = '';
    P = bst_get('ProtocolSubjects');
    if isempty(P)
        error('bst_canonical_cortex:noProtocol', ...
            'No protocol is loaded; cannot resolve the canonical cortex.');
    end
    subj = P.Subject;
    if isfield(P, 'DefaultSubject') && ~isempty(P.DefaultSubject)
        subj = [P.DefaultSubject, subj];
    end

    for k = 1:numel(subj)
        surfs = subj(k).Surface;
        for i = 1:numel(surfs)
            if strcmpi(surfs(i).SurfaceType, 'Cortex')
                m = load(file_fullpath(surfs(i).FileName), 'Vertices');
                if size(m.Vertices, 1) == nVertTarget
                    SurfaceFile = surfs(i).FileName;
                    break;
                end
            end
        end
        if ~isempty(SurfaceFile), break; end
    end

    if isempty(SurfaceFile)
        error('bst_canonical_cortex:notFound', ...
            'No %d-vertex Cortex surface found in the loaded protocol.', nVertTarget);
    end

    % Load mesh data only if the caller asked for it
    if nargout > 1
        TessMat  = in_tess_bst(SurfaceFile);
        Vertices = TessMat.Vertices;
        Faces    = double(TessMat.Faces);
    end
end
