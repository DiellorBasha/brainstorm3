function test_dirac_helmholtz_face()
% Face-NATIVE Helmholtz (bst_helmholtz Domain='face'): the nxr dual gradient gradFace comes
% from the Hodge-Face operator node; geometry (normals/areas/centroids) from the manifold
% node. Convention + exact reconstruction + the strict planted-field round-trip (HarmFrac->0).
% Author: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Surf   = in_tess_bst(SurfaceFile, 0);
    Mani   = tess_manifold(SurfaceFile);
    FaceOp = i_load_op(SurfaceFile, 'Hodge-Face');             % carries FaceAux{hh}.GradFace
    nF = size(Surf.Faces,1);

    Op = bst_helmholtz('Prepare', FaceOp, Mani, Surf, 'Domain','face');
    nFail = nFail + chk('Prepare has G/SkewG + coupled chol per hemi', numel(Op.G)==2 && ...
        size(Op.G{1},1)==3*numel(Op.fH{1}) && size(Op.SkewG{1},2)==numel(Op.fH{1}) && ~isempty(Op.cholA{1}));

    rng(2); Jf = randn(nF, 3) * 1e-9;                          % random per-face field
    Ht = bst_helmholtz('Frame', Op, Jf);

    nFail = nFail + chk('Curl/Div/Psi/Phi are [nF x 1]', isequal(size(Ht.Curl),[nF 1]) && ...
        isequal(size(Ht.Div),[nF 1]) && isequal(size(Ht.Psi),[nF 1]) && isequal(size(Ht.Phi),[nF 1]));
    nFail = nFail + chk('component fields are [nF x 3]', isequal(size(Ht.Virr),[nF 3]) && isequal(size(Ht.Vsol),[nF 3]));
    recon = Ht.Virr + Ht.Vsol + Ht.Vharm;
    nFail = nFail + chk('exact reconstruction Virr+Vsol+Vharm == Jf', max(abs(recon(:)-Ht.Vtot(:))) < 1e-9*max(abs(Ht.Vtot(:))));
    nFail = nFail + chk('HarmFrac in [0,1]', isscalar(Ht.HarmFrac) && Ht.HarmFrac>=0 && Ht.HarmFrac<=1.0001);

    % --- convention: re-decompose the component fields; each must stay in its own channel ---
    en = @(X) sum(X(:).^2);
    Hirr = bst_helmholtz('Frame', Op, Ht.Virr);
    Hsol = bst_helmholtz('Frame', Op, Ht.Vsol);
    nFail = nFail + chk('irrotational stays irrotational (Virr energy > Vsol)', en(Hirr.Virr) > en(Hirr.Vsol));
    nFail = nFail + chk('solenoidal stays solenoidal (Vsol energy > Virr)',     en(Hsol.Vsol) > en(Hsol.Virr));

    nFail = nFail + chk('Cores/Sources are struct arrays w/ persistence', isstruct(Ht.Cores) && isstruct(Ht.Sources) && ...
        (isempty(Ht.Cores) || isfield(Ht.Cores,'persistence')));

    % --- STRICT round-trip: a field planted as SkewG*psi0 must recover exactly ---
    hh=1; fH=Op.fH{hh};
    Cf = Op.Cf{hh};  c = round(size(Cf,1)/2);  d2 = sum((Cf-Cf(c,:)).^2,2);
    psi0 = exp(-d2/(2*0.012^2));  psi0 = psi0 - mean(psi0);
    Vsol0 = reshape(Op.SkewG{hh}*psi0, 3, [])';               % pure solenoidal on faces
    Jsk = zeros(nF,3); Jsk(fH,:) = Vsol0;
    Hsk = bst_helmholtz('Frame', Op, Jsk);
    cc = corr(Hsk.Psi(fH), psi0);
    fprintf('  [planted skew-gradient] HarmFrac=%.3g  corr(psi,psi0)=%.4f\n', Hsk.HarmFrac, cc);
    nFail = nFail + chk('planted skew-gradient: HarmFrac < 0.02', Hsk.HarmFrac < 0.02);
    nFail = nFail + chk('planted skew-gradient: |corr(psi,psi0)| > 0.99', abs(cc) > 0.99);

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
