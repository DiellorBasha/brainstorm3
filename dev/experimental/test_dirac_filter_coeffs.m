function test_dirac_filter_coeffs()
% Projection/filter split equivalence on a tiny synthetic Dirac eigenbasis.
% Authors: Diellor Basha, 2026
    nFail = 0;
    rng(0);
    % two-hemisphere synthetic basis (real eigen nodes always have both hemis):
    % 3 vertices per hemi, K=8 random B-orthonormal Dirac eigenvectors (B = I4 per vertex).
    nVh = 3; K = 8;
    [Q1,~] = qr(randn(4*nVh, K), 0);
    [Q2,~] = qr(randn(4*nVh, K), 0);
    Eigen = struct();
    Eigen.Phi = {Q1, Q2};
    Eigen.Lambda = {sort(rand(K,1)*100), sort(rand(K,1)*100)};
    Eigen.GlobalVertices = {(1:nVh).', (nVh+1:2*nVh).'};
    Mass = {speye(4*nVh), speye(4*nVh)};
    nV = 2*nVh;
    J = randn(3*nV, 1);

    % direct vs split (heat)
    Jdirect = bst_dirac_eigenmodes_filter(Eigen, Mass, J, 'heat', 'DiffusionTime', 0.02);
    [~, ~, c] = bst_dirac_eigenmodes_filter(Eigen, Mass, J, 'custom', 'TransferFn', @(l)ones(size(l)), 'ReturnCoeffs', true);
    Jsplit  = bst_dirac_eigenmodes_filter(Eigen, Mass, [], 'heat', 'DiffusionTime', 0.02, 'Coeffs', c);
    nFail = nFail + chk('split == direct (heat)', max(abs(Jdirect(:)-Jsplit(:))) < 1e-10);

    % re-applying a different gain to the SAME cached coeffs matches a fresh direct call
    Jd2 = bst_dirac_eigenmodes_filter(Eigen, Mass, J, 'bandpass', 'ModeRange', [2 6]);
    Js2 = bst_dirac_eigenmodes_filter(Eigen, Mass, [], 'bandpass', 'ModeRange', [2 6], 'Coeffs', c);
    nFail = nFail + chk('cached coeffs reusable across gains', max(abs(Jd2(:)-Js2(:))) < 1e-10);

    fprintf('\n==== test_dirac_filter_coeffs: %d failed ====\n', nFail);
    if nFail > 0, error('test_dirac_filter_coeffs FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
