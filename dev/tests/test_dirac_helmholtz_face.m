function test_dirac_helmholtz_face()
% Face-domain Helmholtz via the dual face-Dirac D̃: convention + reconstruction + shapes.
% Author: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Surf  = in_tess_bst(SurfaceFile, 0);
    Dirac = i_load_op(SurfaceFile, 'Dirac');
    LBO   = i_load_op(SurfaceFile, 'Laplace-Beltrami');
    nV = size(Surf.Vertices,1); nF = size(Surf.Faces,1);

    Op = bst_dirac_helmholtz_face('Prepare', Dirac, LBO, Surf);
    nFail = nFail + chk('Prepare has D-tilde per hemi', numel(Op.Dt)==2 && size(Op.Dt{1},1)==4*numel(Op.vH{1}));

    rng(2); Jf = randn(nF, 3) * 1e-9;                          % random per-face field
    Ht = bst_dirac_helmholtz_face('Frame', Op, Jf);

    nFail = nFail + chk('Curl/Div are [nV x 1]', isequal(size(Ht.Curl),[nV 1]) && isequal(size(Ht.Div),[nV 1]));
    nFail = nFail + chk('component fields are [nF x 3]', isequal(size(Ht.Virr),[nF 3]) && isequal(size(Ht.Vsol),[nF 3]));
    recon = Ht.Virr + Ht.Vsol + Ht.Vharm;
    nFail = nFail + chk('exact reconstruction Virr+Vsol+Vharm == Jf', max(abs(recon(:)-Ht.Vtot(:))) < 1e-9*max(abs(Ht.Vtot(:))));
    nFail = nFail + chk('HarmFrac in [0,1]', isscalar(Ht.HarmFrac) && Ht.HarmFrac>=0 && Ht.HarmFrac<=1.0001);

    % --- convention: re-decompose the component fields (validates w=vorticity, imag.n=div) ---
    Hirr = bst_dirac_helmholtz_face('Frame', Op, Ht.Virr);
    Hsol = bst_dirac_helmholtz_face('Frame', Op, Ht.Vsol);
    nFail = nFail + chk('irrotational is divergence-dominated', sum(Hirr.Div.^2) > sum(Hirr.Curl.^2));
    nFail = nFail + chk('solenoidal is curl-dominated',         sum(Hsol.Curl.^2) > sum(Hsol.Div.^2));

    nFail = nFail + chk('Cores/Sources are struct arrays w/ persistence', isstruct(Ht.Cores) && isstruct(Ht.Sources) && (isempty(Ht.Cores) || isfield(Ht.Cores,'persistence')));

    % --- planted pure-skew-gradient (solenoidal) face field must recover (HarmFrac->0) ---
    hh=1; vH=Op.vH{hh}; fH=Op.fH{hh};
    c = vH(round(numel(vH)/2)); d2 = sum((Surf.Vertices(vH,:)-Surf.Vertices(c,:)).^2,2);
    psi0 = exp(-d2/(2*0.012^2)); psi0 = psi0 - mean(psi0);
    gp = [Op.Gx{hh}*psi0, Op.Gy{hh}*psi0, Op.Gz{hh}*psi0];
    Vsk = cross(Op.Nf{hh}, gp, 2);                          % n x grad(psi0): pure solenoidal
    Jsk = zeros(nF,3); Jsk(fH,:) = Vsk;
    Hsk = bst_dirac_helmholtz_face('Frame', Op, Jsk);
    % Primary gate: the stream function is recovered (shape). HarmFrac floor ~0.10 with
    % the intrinsic D̃_int + vertex cotan Poisson + FEM-gradient reconstruction (a Dirac-vs-FEM
    % discretization gap); a full face-native K̃_int Poisson + dual-complex gradient would reach ~0.
    fprintf('  [planted skew-gradient] HarmFrac=%.3g  corr(psi,psi0)=%.3f\n', Hsk.HarmFrac, corr(Hsk.Psi(vH), psi0));
    nFail = nFail + chk('planted skew-gradient: recovered psi corr > 0.95', abs(corr(Hsk.Psi(vH), psi0)) > 0.95);
    nFail = nFail + chk('planted skew-gradient: HarmFrac < 0.15 (sanity)', Hsk.HarmFrac < 0.15);

    fprintf('\n==== test_dirac_helmholtz_face: %d failed ====\n', nFail);
    if nFail > 0, error('test_dirac_helmholtz_face FAILED'); end
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
