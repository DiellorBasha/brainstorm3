function test_dirac_face_operator_bst()
% Phase 1: tess_operators 'Dirac-Face' variant -- A=(1-Tau)E~_int+Tau.E~_ext [4F x 4F],
% mass W_F, co-normalized. Pins the W_V metric against nxr extrinsicBlockFace.
% Author: Diellor Basha, 2026
    nFail = 0;  Tau = 0.5;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Op = tess_operators(SurfaceFile, 'Dirac-Face', 'Tau', Tau, 'NoSave', true);
    nFail = nFail + chk('Variant tag', strcmpi(Op.Variant,'Dirac-Face'));
    for hh = 1:2
        A = Op.Operator{hh}; B = Op.Mass{hh}; n4 = size(A,1);
        nFail = nFail + chk(sprintf('h%d A,B square & equal size',hh), size(A,1)==size(A,2) && isequal(size(B),size(A)));
        nFail = nFail + chk(sprintf('h%d A symmetric',hh), norm(A-A','fro') < 1e-8*norm(A,'fro'));
        % A is a positive combination of two Galerkin squares (D~'W D~) -> PSD by construction.
        % Verify via the min generalized Rayleigh quotient over random probes (bulletproof; no
        % eigensolver convergence dependence -- the full spectrum is validated in Phase 2).
        rng(7); X = randn(size(A,1), 40);
        rq = sum(X.*(A*X),1) ./ max(sum(X.*(B*X),1), eps);
        fprintf('  h%d min Rayleigh quotient over 40 probes = %.3e\n', hh, min(rq));
        nFail = nFail + chk(sprintf('h%d A PSD (min Rayleigh quotient >= -1e-9)',hh), min(rq) >= -1e-9);
        % the smallest TRUE eigenvalues (shift-invert at a small positive sigma) are ~0 (constants)
        d = sort(real(eigs(0.5*(A+A'), 0.5*(B+B'), 8, 1e-6)));
        fprintf('  h%d smallest 8 eigs: %s\n', hh, mat2str(d',3));
        nFail = nFail + chk(sprintf('h%d smallest eig >= -1e-6 (PSD spectrum)',hh), d(1) >= -1e-6);
        nFail = nFail + chk(sprintf('h%d B diagonal (W_F)',hh), nnz(B - diag(diag(B)))==0);
        nFail = nFail + chk(sprintf('h%d GlobalFaces 4*nF == size(A)',hh), ~isempty(Op.GlobalFaces{hh}) && 4*numel(Op.GlobalFaces{hh})==n4);
    end

    % --- W_V metric pin: assembling E~_ext ourselves with W_V (vertex dual area, lumped
    %     mass) must reproduce nxr extrinsicBlockFace. This pins the metric E~_int shares. ---
    [Vh,Fh] = i_hemi_submesh(SurfaceFile, 1);
    h = nxr_safe_create(Vh, Fh);
    Dext     = nxr_compute('operators', h, 'diracFaceD');     % [4V x 4F]
    Eext_nxr = nxr_compute('operators', h, 'diracFace', 1);   % extrinsicBlockFace [4F x 4F]
    Ml       = nxr_compute('operators', h, 'mass', 'lumped'); % vertex dual area [nV x nV]
    nxr_compute('destroy', h);
    WV = kron(Ml, speye(4));
    Eext_mine = Dext' * WV * Dext;
    rel = norm(Eext_mine - Eext_nxr,'fro') / max(norm(Eext_nxr,'fro'),eps);
    fprintf('  [W_V metric pin] ||D~_ext'' W_V D~_ext - extrinsicBlockFace|| / ||.|| = %.2e\n', rel);
    nFail = nFail + chk('W_V (lumped) metric pin < 1e-7', rel < 1e-7);

    fprintf('\n==== test_dirac_face_operator_bst: %d failed ====\n', nFail);
    if nFail > 0, error('test_dirac_face_operator_bst FAILED'); end
end

function [Vloc, Floc] = i_hemi_submesh(SurfaceFile, hh)
% Build hemisphere-hh local submesh exactly as tess_operators does (tess_hemisplit).
    TessMat = in_tess_bst(SurfaceFile, 0);
    [ir, il] = tess_hemisplit(TessMat);
    hemis = {il, ir};  vH = double(hemis{hh}(:));
    Vtx = TessMat.Vertices; Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);
    isV = false(nVtot,1); isV(vH)=true; fMask = all(isV(Fcs),2);
    mapV = zeros(nVtot,1); mapV(vH)=1:numel(vH);
    Vloc = Vtx(vH,:); Floc = mapV(Fcs(fMask,:));
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
