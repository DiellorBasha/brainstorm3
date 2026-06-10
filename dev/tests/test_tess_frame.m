function tests = test_tess_frame
% Tests for tess_frame: compute/store the nxr facet bundle + return {U,V,N}.
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

function test_store_five_groups_1x2(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for f = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'}
        verifyTrue(tc, isfield(T, f{1}) && isequal(size(T.(f{1})), [1 2]), ...
            sprintf('%s is 1x2', f{1}));
    end
    verifyEqual(tc, T.Embedded(1).Hemisphere, 'L');
    verifyEqual(tc, T.Embedded(2).Hemisphere, 'R');
    verifyEqual(tc, T.Gauge(1).Provenance.Backend, 'nxr');
    verifyEqual(tc, lower(T.Gauge(1).type), 'trivial');
end

function test_scatter_partition(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    gv = sort([T.Embedded(1).GlobalVertices(:); T.Embedded(2).GlobalVertices(:)]);
    gf = sort([T.Embedded(1).GlobalFaces(:);    T.Embedded(2).GlobalFaces(:)]);
    verifyEqual(tc, gv, (1:size(T.Vertices,1))');
    verifyEqual(tc, gf, (1:size(T.Faces,1))');
end

function test_frame_orthonormal_fullmesh(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    [U,V,N] = tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    nV = size(T.Vertices,1);
    verifyEqual(tc, size(U), [nV 3]);
    verifyEqual(tc, size(V), [nV 3]);
    verifyEqual(tc, size(N), [nV 3]);
    verifyLessThan(tc, max(abs(sqrt(sum(U.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sqrt(sum(V.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sum(U.*V,2))), 1e-4);
    verifyLessThan(tc, max(abs(N - cross(U,V,2)), [], 'all'), 1e-4);
    verifyGreaterThan(tc, min(sqrt(sum(N.^2,2))), 0.9);   % every row filled
end

function test_frame_matches_grid_rotation(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    [U,V,~] = tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        idx = T.Embedded(hh).GlobalVertices;
        c   = T.Embedded(hh).vertex.grid .* T.Gauge(hh).vertex.rotation;
        verifyLessThan(tc, max(abs(U(idx,:) - real(c)), [], 'all'), 1e-5);
        verifyLessThan(tc, max(abs(V(idx,:) - imag(c)), [], 'all'), 1e-5);
    end
end
