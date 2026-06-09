function tests = test_tess_topology
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

function test_topology_struct_array(tc)
    Topo = tess_topology(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(Topo), [1 2]);
    verifyTrue(tc, isfield(Topo(1), 'schemaVersion'));
    verifyEqual(tc, Topo(1).Hemisphere, 'L');
    verifyEqual(tc, Topo(2).Hemisphere, 'R');
    verifyEqual(tc, Topo(1).Provenance.Package, 'topology');
end

function test_topology_partition(tc)
    SurfaceFile = local_cortex();
    Topo = tess_topology(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    gv = sort([Topo(1).GlobalVertices(:); Topo(2).GlobalVertices(:)]);
    verifyEqual(tc, gv, (1:size(T.Vertices,1))');
end

function test_topology_operators(tc)
    Topo = tess_topology(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1);
    for hh = 1:2
        op = Topo(hh).operators;
        verifyEqual(tc, nnz(op.dec.d1 * op.dec.d0), 0);   % d o d = 0 (exact)
    end
end

function test_topology_is_stored(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_topology(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Topology') && isequal(size(T.Topology),[1 2]));
end
