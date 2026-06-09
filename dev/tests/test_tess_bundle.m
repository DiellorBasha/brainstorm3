function tests = test_tess_bundle
% Property tests for tess_bundle (per-hemisphere nxr bundle store) on real cortex.
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

function test_three_struct_arrays_1x2(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(B.Topology), [1 2]);
    verifyEqual(tc, size(B.Geometry), [1 2]);
    verifyEqual(tc, size(B.Gauge),    [1 2]);
    verifyEqual(tc, B.Gauge(1).Provenance.Backend, 'nxr');   % writer-added provenance
end

function test_scatter_maps_partition_mesh(tc)
    SurfaceFile = local_cortex();
    B = tess_bundle(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    gv = [B.Geometry(1).GlobalVertices(:); B.Geometry(2).GlobalVertices(:)];
    gf = [B.Geometry(1).GlobalFaces(:);    B.Geometry(2).GlobalFaces(:)];
    % disjoint + complete partition of vertices and faces
    verifyEqual(tc, sort(gv), (1:size(T.Vertices,1))');
    verifyEqual(tc, sort(gf), (1:size(T.Faces,1))');
    verifyEqual(tc, B.Geometry(1).Hemisphere, 'L');
    verifyEqual(tc, B.Geometry(2).Hemisphere, 'R');
end

function test_no_operators_in_light(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyFalse(tc, isfield(B.Topology(1), 'operators') && ~isempty(B.Topology(1).operators));
    verifyFalse(tc, isfield(B.Geometry(1), 'operators') && ~isempty(B.Geometry(1).operators));
    verifyFalse(tc, isfield(B.Gauge(1),    'operators') && ~isempty(B.Gauge(1).operators));
end

function test_vertex_grid_orthonormal(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    for hh = 1:2
        c  = B.Geometry(hh).vertex.grid;        % nVh x 3 complex
        e1 = real(c); e2 = imag(c);
        n1 = sqrt(sum(e1.^2,2)); n2 = sqrt(sum(e2.^2,2));
        verifyLessThan(tc, max(abs(n1-1)), 1e-4);
        verifyLessThan(tc, max(abs(n2-1)), 1e-4);
        verifyLessThan(tc, max(abs(sum(e1.*e2,2))), 1e-4);
    end
end

function test_trivial_gauge_gauss_bonnet(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);   % default gauge='trivial'
    for hh = 1:2
        verifyEqual(tc, lower(B.Gauge(hh).type), 'trivial');
        verifyEqual(tc, sum(B.Gauge(hh).singularity.indices), 2, 'AbsTol', 1e-6);
    end
end

function test_field_is_stored(tc)
    % Real save path, then restore the file (isolation).
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_bundle(SurfaceFile, 'ForceRecompute', 1);   % default: save
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Topology') && isequal(size(T.Topology),[1 2]));
    verifyTrue(tc, isfield(T,'Geometry') && isequal(size(T.Geometry),[1 2]));
    verifyTrue(tc, isfield(T,'Gauge')    && isequal(size(T.Gauge),[1 2]));
end
