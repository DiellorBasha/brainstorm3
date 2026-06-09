function tests = test_tess_geometry
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

function test_geometry_struct_array(tc)
    Geo = tess_geometry(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(Geo), [1 2]);
    verifyEqual(tc, Geo(1).Provenance.Package, 'geometry');
end

function test_geometry_grid_orthonormal(tc)
    Geo = tess_geometry(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    for hh = 1:2
        c = Geo(hh).vertex.grid; e1 = real(c); e2 = imag(c);
        verifyLessThan(tc, max(abs(sqrt(sum(e1.^2,2))-1)), 1e-4);
        verifyLessThan(tc, max(abs(sqrt(sum(e2.^2,2))-1)), 1e-4);
        verifyLessThan(tc, max(abs(sum(e1.*e2,2))), 1e-4);
    end
end

function test_geometry_operators(tc)
    Geo = tess_geometry(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1, 'Mass', 'lumped');
    for hh = 1:2
        op = Geo(hh).operators;
        verifyTrue(tc, isfield(op,'laplacian') && isfield(op,'mass') && isfield(op,'hodge'));
        verifyTrue(tc, issparse(op.mass.lumped));
    end
end

function test_geometry_is_stored(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_geometry(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Geometry') && isequal(size(T.Geometry),[1 2]));
end
