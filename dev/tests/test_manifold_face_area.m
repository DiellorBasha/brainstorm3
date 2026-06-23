function test_manifold_face_area()
% TEST_MANIFOLD_FACE_AREA: the manifold node carries Embedded.face.area (nxr facets
% schemaVersion >= 2), and tess_manifold's schema gate recomputes a stale (pre-area)
% node instead of returning it.
%
% Needs a loaded protocol whose Subject 1 has a cortex surface; SKIPs otherwise.
% Author: Diellor Basha, 2026
    if isempty(bst_get('ProtocolInfo'))
        disp('SKIP: test_manifold_face_area -- no protocol loaded'); return;
    end
    sSubject = bst_get('Subject', 1);
    if isempty(sSubject) || isempty(sSubject.iCortex)
        disp('SKIP: test_manifold_face_area -- Subject 1 has no cortex'); return;
    end
    SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;

    % --- 1) compute-and-save a fresh (schemaVersion 2) manifold ---
    M = tess_manifold(SurfaceFile);                 % find-or-create + save
    assert(M.Embedded(1).schemaVersion >= 2 && isfield(M.Embedded(1).face,'area'), ...
        'fresh manifold lacks schemaVersion>=2 / face.area');
    for hh = 1:numel(M.Embedded)
        a = M.Embedded(hh).face.area;
        assert(numel(a)==numel(M.Embedded(hh).GlobalFaces) && all(a>0), 'h%d area invalid', hh);
        fprintf('  fresh h%d nF=%d sum(area)=%.5f m^2\n', hh, numel(a), sum(a));
    end

    % --- 2) downgrade the on-disk node to simulate a pre-area (schema 1) node ---
    [sMan, ~, iSurfMan, iMan] = bst_get('ManifoldFileForSurface', SurfaceFile, 'trivial');
    assert(~isempty(iMan), 'no manifold node registered after save');
    existFile = file_fullpath(sMan.Surface(iSurfMan).Manifold(iMan).FileName);
    D = load(existFile);
    for hh = 1:numel(D.Embedded)
        D.Embedded(hh).schemaVersion = 1;
        if isfield(D.Embedded(hh).face,'area')
            D.Embedded(hh).face = rmfield(D.Embedded(hh).face, 'area');
        end
    end
    bst_save(existFile, D, 'v7');
    % confirm the on-disk node is now stale
    Dchk = in_bst_manifold(existFile);
    assert(Dchk.Embedded(1).schemaVersion==1 && ~isfield(Dchk.Embedded(1).face,'area'), ...
        'downgrade did not take');
    disp('  on-disk node downgraded to schemaVersion 1 (no face.area)');

    % --- 3) tess_manifold WITHOUT ForceRecompute must auto-recompute the stale node ---
    M2 = tess_manifold(SurfaceFile);
    assert(M2.Embedded(1).schemaVersion >= 2 && isfield(M2.Embedded(1).face,'area'), ...
        'GATE FAILED: stale node was returned instead of recomputed');
    for hh = 1:numel(M2.Embedded)
        assert(numel(M2.Embedded(hh).face.area)==numel(M2.Embedded(hh).GlobalFaces) ...
            && all(M2.Embedded(hh).face.area>0), 'h%d recomputed area invalid', hh);
    end
    disp('PASS: test_manifold_face_area (face.area present; stale-schema node auto-recomputed)');
end
