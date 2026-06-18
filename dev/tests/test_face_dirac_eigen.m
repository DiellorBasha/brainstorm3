function test_face_dirac_eigen()
% Phase 2: tess_eigen 'Dirac-Face' eigenmodes Phi [4nF x K] via the SAME machinery
% as the vertex Dirac (bst_eigs_smallest + Rayleigh-Ritz). Gate = B-orthonormality.
% Author: Diellor Basha, 2026
    nFail = 0;  K = 200;  Tau = 0.5;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    Eig = tess_eigen(SurfaceFile, 'Dirac-Face', 'K', K, 'Tau', Tau, 'NoSave', true);
    Op  = tess_operators(SurfaceFile, 'Dirac-Face', 'Tau', Tau, 'NoSave', true);
    nFail = nFail + chk('Variant tag', strcmpi(Eig.Variant,'Dirac-Face'));
    for hh = 1:2
        Phi = Eig.Phi{hh};  B = Op.Mass{hh};  lam = Eig.Lambda{hh}(:);
        nFail = nFail + chk(sprintf('h%d Phi [4nF x K]',hh), size(Phi,1)==size(B,1) && size(Phi,2)>=K);
        G = Phi(:,1:K)' * B * Phi(:,1:K);
        ortho = norm(G - eye(K), 'fro') / sqrt(K);
        fprintf('  h%d B-orthonormality ||Phi''BPhi - I||/sqrt(K) = %.2e\n', hh, ortho);
        nFail = nFail + chk(sprintf('h%d B-orthonormal',hh), ortho < 1e-5);
        fprintf('  h%d Lambda(1:6) = %s ... Lambda(K)=%.3g\n', hh, mat2str(lam(1:6)',3), lam(min(K,end)));
        nFail = nFail + chk(sprintf('h%d Lambda ascending, >= -1e-6',hh), all(diff(lam(1:K))>=-1e-8) && min(lam(1:K))>=-1e-6);
        nFail = nFail + chk(sprintf('h%d GlobalFaces set & matches 4nF',hh), ...
            ~isempty(Eig.GlobalFaces{hh}) && 4*numel(Eig.GlobalFaces{hh})==size(Phi,1));
    end
    fprintf('\n==== test_face_dirac_eigen: %d failed ====\n', nFail);
    if nFail > 0, error('test_face_dirac_eigen FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
