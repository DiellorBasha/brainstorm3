function test_bst_eigen_resolve()
% TEST_BST_EIGEN_RESOLVE: bst_eigen resolves the eigen_ basis implicitly from a results
% file's SurfaceFile (Variant=Laplace-Beltrami), and raises a clear error when none exists.
% SKIPs if the live protocol lacks an LBO eigen node or a constrained results file on its surface.
    nFail = 0; chk = @i_chk;
    [ef, surf, nv] = i_find_eigen('Laplace-Beltrami');
    if isempty(ef); fprintf('SKIP test_bst_eigen_resolve: no LBO eigen node.\n'); return; end
    fprintf('LBO eigen: %s\n  surface: %s\n', ef, surf);

    % (1) resolver returns a node for the surface
    [sS, ~, iSurf, iE] = bst_get('EigenFileForSurface', surf, 'Laplace-Beltrami');
    nFail = nFail + chk('EigenFileForSurface resolves LBO', ~isempty(iE) && ~isempty(sS));
    if ~isempty(iE)
        nFail = nFail + chk('resolved entry has a FileName', ~isempty(sS.Surface(iSurf).Eigen(iE).FileName));
    end

    % (2) end-to-end implicit resolution on a real constrained results file
    rf = i_find_results(surf, nv);
    if isempty(rf)
        fprintf('SKIP end-to-end: no constrained results on surface.\n');
    else
        fprintf('results: %s\n', rf);
        O = bst_eigen(); O.Method = 'spectrum'; O.Variant = 'Laplace-Beltrami';
        O.EigenFile = []; O.iTargetStudy = 'NoSave';
        [out, msg, err] = bst_eigen(rf, O);
        nFail = nFail + chk('implicit resolve: no error', err == 0 && isempty(msg));
        ok = iscell(out) && ~isempty(out) && isstruct(out{1});
        nFail = nFail + chk('implicit resolve: returns a struct', ok);
        if ok
            nFail = nFail + chk('spectrum has Freqs (sqrt-lambda)', isfield(out{1}, 'Freqs') && ~isempty(out{1}.Freqs));
            nFail = nFail + chk('spectrum Method == spectrum', isfield(out{1}, 'Method') && strcmp(out{1}.Method, 'spectrum'));
        end
    end

    % (3) negative: a variant absent on this surface -> clear typed error
    bad = 'Connection Laplacian';
    [~, ~, ~, iEb] = bst_get('EigenFileForSurface', surf, bad);
    if isempty(iEb) && ~isempty(rf)
        Ob = bst_eigen(); Ob.Method = 'spectrum'; Ob.Variant = bad;
        Ob.EigenFile = []; Ob.iTargetStudy = 'NoSave';
        gotErr = false;
        try
            bst_eigen(rf, Ob);
        catch ME
            gotErr = strcmp(ME.identifier, 'bst_eigen:NoEigenForSurface');
        end
        nFail = nFail + chk('missing variant -> NoEigenForSurface error', gotErr);
    else
        fprintf('SKIP negative: a Connection basis exists or no results.\n');
    end

    fprintf('\n==== test_bst_eigen_resolve: %d failed ====\n', nFail);
    if nFail > 0; error('test_bst_eigen_resolve FAILED'); end
    disp('ALL TESTS PASSED');
end

function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end
function [ef, surf, nv] = i_find_eigen(want)
    % Returns the most recent eigen node of the wanted Variant, its parent surface, and the
    % basis vertex count nv (max GlobalVertices) used to dimension-match a compatible results file.
    ef = ''; surf = ''; nv = 0;
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat')); [~, o] = sort([d.datenum], 'descend');
    for i = o(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant', 'ParentSurface', 'GlobalVertices');
            if strcmpi(m.Variant, want)
                ef = rel; surf = m.ParentSurface;
                for h = 1:numel(m.GlobalVertices)
                    if ~isempty(m.GlobalVertices{h}); nv = max(nv, max(m.GlobalVertices{h})); end
                end
                return;
            end
        catch
        end
    end
end
function rf = i_find_results(surf, expectNV)
    % First constrained (nComponents==1) results file on surf whose source count == expectNV
    % (the basis vertex count). Surface NAME alone is insufficient: a stale results file from
    % before the surface was re-tessellated can share the name but have a different vertex count.
    rf = '';
    sP = bst_get('ProtocolStudies'); if isempty(sP); return; end
    for s = 1:numel(sP.Study)
        R = sP.Study(s).Result;
        for r = 1:numel(R)
            try
                m = in_bst_results(R(r).FileName, 0, 'SurfaceFile', 'nComponents');
                if isempty(m.SurfaceFile) || ~file_compare(m.SurfaceFile, surf) ...
                        || ~isequal(m.nComponents, 1)
                    continue;
                end
                mf = in_bst_results(R(r).FileName, 1, 'ImageGridAmp');
                if size(mf.ImageGridAmp, 1) == expectNV
                    rf = R(r).FileName; return;
                end
            catch
            end
        end
    end
end
