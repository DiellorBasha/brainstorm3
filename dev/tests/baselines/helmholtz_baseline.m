function helmholtz_baseline
% Capture old bst_helmholtz('Frame') outputs as the refactor parity oracle.
% Run ONCE on the unmodified tree, BEFORE editing any differential/ file.
% HISTORICAL: bst_helmholtz was deleted at the end of the refactor, so this
% capture script no longer runs. The oracle it produced (helmholtz_baseline.mat)
% is committed and is what the parity tests load. Kept for provenance only.
    SurfaceFile = 'Subject01/tess_cortex_pial_low.mat';   % cortex w/ manifold_ + Covariant ready
    Surf = in_tess_bst(SurfaceFile, 0);
    nV   = size(Surf.Vertices, 1);
    Cov  = tess_operators(SurfaceFile, 'Covariant');
    LBO  = tess_operators(SurfaceFile, 'Laplace-Beltrami');
    Mani = tess_manifold(SurfaceFile);
    rng(7);  J = randn(3*nV, 1);                          % fixed-seed synthetic source frame

    Op = bst_helmholtz('Prepare', {Cov, LBO}, Mani, Surf, 'Domain', 'vertex');
    Ht = bst_helmholtz('Frame', Op, J, false);            % cores off

    B = struct('Surf', SurfaceFile, 'J', J, ...
        'Div', Ht.Div, 'Curl', Ht.Curl, 'Phi', Ht.Phi, 'Psi', Ht.Psi, ...
        'Fmag', Ht.Fmag, 'Hmag', Ht.Hmag, 'Virr', Ht.Virr, 'Vsol', Ht.Vsol, ...
        'Vtot', Ht.Vtot, 'Hresid', Ht.Hresid, 'HarmFrac', Ht.HarmFrac); %#ok<NASGU>
    outDir = fileparts(mfilename('fullpath'));
    save(fullfile(outDir, 'helmholtz_baseline.mat'), 'B');
    fprintf('helmholtz_baseline: saved %d-vertex baseline (HarmFrac=%.3e)\n', nV, Ht.HarmFrac);
end
