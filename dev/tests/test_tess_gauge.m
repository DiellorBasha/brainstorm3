function tests = test_tess_gauge
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

function test_gauge_struct_array(tc)
    Ga = tess_gauge(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);   % default trivial
    verifyEqual(tc, size(Ga), [1 2]);
    verifyEqual(tc, lower(Ga(1).type), 'trivial');
    verifyEqual(tc, Ga(1).Provenance.Package, 'gauge');
end

function test_gauge_gauss_bonnet(tc)
    Ga = tess_gauge(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    for hh = 1:2
        verifyEqual(tc, sum(Ga(hh).singularity.indices), 2, 'AbsTol', 1e-6);
    end
end

function test_gauge_operators(tc)
    Ga = tess_gauge(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1, 'Coupling', 'ambient');
    for hh = 1:2
        nVh = numel(Ga(hh).GlobalVertices);
        K = Ga(hh).operators.laplacian;
        verifyTrue(tc, issparse(K));
        verifyGreaterThan(tc, full(max(abs(imag(K)), [], 'all')), 0);     % genuinely complex
        verifyLessThan(tc, full(max(abs(K - K'), [], 'all')), 1e-6);      % Hermitian
        L3 = Ga(hh).operators.covariantLaplacian;
        verifyEqual(tc, size(L3), [3*nVh, 3*nVh]);
    end
end

function test_gauge_is_stored(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_gauge(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Gauge') && isequal(size(T.Gauge),[1 2]));
end
