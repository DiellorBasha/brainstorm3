function tests = test_tess_frame
% Tests for tess_frame (derived {U,V,N} from the stored bundle) on real cortex.
tests = functiontests(localfunctions);
end

function SurfaceFile = local_cortex()
    if ~brainstorm('status'); brainstorm nogui; end
    SurfaceFile = local_find_cortex(20484);
end

function s = local_find_cortex(nVertTarget)
    s = '';
    P = bst_get('ProtocolSubjects');
    subj = P.Subject;
    if isfield(P,'DefaultSubject') && ~isempty(P.DefaultSubject)
        subj = [P.DefaultSubject, subj];
    end
    for k = 1:numel(subj)
        surfs = subj(k).Surface;
        for i = 1:numel(surfs)
            if strcmpi(surfs(i).SurfaceType, 'Cortex')
                m = load(file_fullpath(surfs(i).FileName), 'Vertices');
                if size(m.Vertices,1) == nVertTarget
                    s = surfs(i).FileName; return;
                end
            end
        end
    end
    error('No %d-vertex cortex found in the loaded protocol.', nVertTarget);
end

function local_ensure_bundle(SurfaceFile)
    % Make sure a trivial-gauge bundle is on the file; restore is the caller's job.
    T = in_tess_bst(SurfaceFile, 0);
    if ~isfield(T,'Geometry') || isempty(T.Geometry)
        tess_bundle(SurfaceFile, 'ForceRecompute', 1);
    end
end

function test_frame_fullmesh_orthonormal(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    local_ensure_bundle(SurfaceFile);

    [U,V,N] = tess_frame(SurfaceFile);          % default vertex domain
    T = in_tess_bst(SurfaceFile, 0);
    nV = size(T.Vertices,1);
    verifyEqual(tc, size(U), [nV 3]);
    verifyEqual(tc, size(V), [nV 3]);
    verifyEqual(tc, size(N), [nV 3]);
    % orthonormal, right-handed, no unfilled rows
    verifyLessThan(tc, max(abs(sqrt(sum(U.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sqrt(sum(V.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sum(U.*V,2))), 1e-4);
    verifyLessThan(tc, max(abs(N - cross(U,V,2)), [], 'all'), 1e-4);
    verifyGreaterThan(tc, min(sqrt(sum(N.^2,2))), 0.9);    % every row filled
end

function test_frame_matches_grid_rotation(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    local_ensure_bundle(SurfaceFile);

    [U,V,~] = tess_frame(SurfaceFile);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        idx = T.Geometry(hh).GlobalVertices;
        c   = T.Geometry(hh).vertex.grid .* T.Gauge(hh).vertex.rotation;   % grid (x) rotation
        verifyLessThan(tc, max(abs(U(idx,:) - real(c)), [], 'all'), 1e-5);
        verifyLessThan(tc, max(abs(V(idx,:) - imag(c)), [], 'all'), 1e-5);
    end
end

function test_face_trivial_errors(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    local_ensure_bundle(SurfaceFile);
    verifyError(tc, @() tess_frame(SurfaceFile, 'Domain','face'), 'tess_frame:faceTrivialDeferred');
end
