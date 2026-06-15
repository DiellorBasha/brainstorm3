function test_dirac_helmholtz()
% Pure tests: Poisson solve round-trip, core classifier on a planted scalar, and a
% full-pipeline sanity on the real Subject01 cortex (Dirac + LBO operators).
% Authors: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Surf  = in_tess_bst(SurfaceFile, 0);
    Dirac = i_load_op(SurfaceFile, 'Dirac');
    LBO   = i_load_op(SurfaceFile, 'Laplace-Beltrami');

    % --- (1) Poisson solve round-trip on hemisphere 1: K psi = M omega recovers psi ---
    K = LBO.Operator{1}; M = LBO.Mass{1}; n = size(K,1);
    rng(0); psiTrue = randn(n,1); psiTrue = psiTrue - mean(psiTrue);
    rhs = K * psiTrue;                                   % = M*omega with omega = Lap*psiTrue
    psiRec = bst_dirac_helmholtz('PoissonSolve', K, M, M \ rhs);  % omega = M\rhs
    psiRec = psiRec - mean(psiRec);
    nFail = nFail + chk('Poisson recovers psi (up to const)', max(abs(psiRec - psiTrue)) < 1e-6 * max(abs(psiTrue)));

    % --- (2) Core classifier: a Gaussian bump -> its peak is detected as a (max) core ---
    % (A single bump necessarily also yields a min/saddle elsewhere -- the index theorem;
    %  the essential property is that the planted peak is found and classified as a maximum.)
    V = Surf.Vertices; vc = 5000;
    d2 = sum((V - V(vc,:)).^2, 2);
    psi = exp(-d2 / (2*(0.01)^2));                      % sharp positive bump at vc
    cores = bst_dirac_helmholtz('FindCores', psi, Surf.VertConn, zeros(size(psi)));
    iVs = [cores.iVertex];
    nFail = nFail + chk('bump peak is a core', ismember(vc, iVs));
    cvc = cores(iVs==vc);
    nFail = nFail + chk('peak classified as a maximum (charge<0 w/ omega=0)', ~isempty(cvc) && cvc(1).charge < 0);

    % --- (3) Full pipeline on a real (random) source field: finite, right size, nonzero ---
    nV = size(V,1); rng(1); J = randn(3*nV, 2) * 1e-9;
    H = bst_dirac_helmholtz(Dirac, LBO, Surf, J);
    nFail = nFail + chk('Curl size [nV x nT]', isequal(size(H.Curl), [nV 2]));
    nFail = nFail + chk('Div/Psi/Phi present + finite', all(isfinite(H.Div(:))) && all(isfinite(H.Psi(:))) && all(isfinite(H.Phi(:))));
    nFail = nFail + chk('field has both curl and div', max(abs(H.Curl(:)))>0 && max(abs(H.Div(:)))>0);
    nFail = nFail + chk('Cores is 1xnT cell', iscell(H.Cores) && numel(H.Cores)==2);

    % --- (4) on-demand path: Prepare once + Frame(single col) == whole-series column ---
    Op  = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    Ht  = bst_dirac_helmholtz('Frame', Op, J(:,1));
    nFail = nFail + chk('Frame Curl == Decompose col 1', isequal(Ht.Curl, H.Curl(:,1)));
    nFail = nFail + chk('Frame Psi  == Decompose col 1', isequal(Ht.Psi,  H.Psi(:,1)));
    nFail = nFail + chk('Frame cores match col 1', numel(Ht.Cores)==numel(H.Cores{1}));

    % --- (5) Hodge decomposition: component fields, exact reconstruction, dominance ---
    Op2 = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    Ht  = bst_dirac_helmholtz('Frame', Op2, J(:,1));
    nFail = nFail + chk('component fields are [nV x 3]', isequal(size(Ht.Virr),[size(V,1) 3]) && isequal(size(Ht.Vsol),[size(V,1) 3]));
    recon = Ht.Virr + Ht.Vsol + Ht.Vharm;
    nFail = nFail + chk('exact reconstruction Virr+Vsol+Vharm == J', max(abs(recon(:)-Ht.Vtot(:))) < 1e-9*max(abs(Ht.Vtot(:))));
    Hirr = bst_dirac_helmholtz('Frame', Op2, reshape(Ht.Virr',[],1));
    Hsol = bst_dirac_helmholtz('Frame', Op2, reshape(Ht.Vsol',[],1));
    nFail = nFail + chk('irrotational is divergence-dominated', sum(Hirr.Div.^2) > sum(Hirr.Curl.^2));
    nFail = nFail + chk('solenoidal is curl-dominated',         sum(Hsol.Curl.^2) > sum(Hsol.Div.^2));
    nFail = nFail + chk('Cores + Sources are struct arrays', isstruct(Ht.Cores) && isstruct(Ht.Sources));
    nFail = nFail + chk('HarmFrac in [0,1]', isscalar(Ht.HarmFrac) && Ht.HarmFrac >= 0 && Ht.HarmFrac <= 1.0001);

    fprintf('\n==== test_dirac_helmholtz: %d failed ====\n', nFail);
    if nFail > 0, error('test_dirac_helmholtz FAILED'); end
end

function Op = i_load_op(SurfaceFile, variant)
    sSubject = bst_get('Subject',1);
    iSurf = find(strcmpi({sSubject.Surface.FileName}, file_short(SurfaceFile)),1);
    Op = [];
    if isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = load(file_fullpath(sSubject.Surface(iSurf).Operator(k).FileName));
            if strcmpi(S.Variant, variant); Op = S; break; end
        end
    end
    if isempty(Op); tess_operators(SurfaceFile, variant); Op = i_load_op(SurfaceFile, variant); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
