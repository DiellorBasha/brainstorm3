function test_process_eigen()
% TEST_PROCESS_EIGEN: the pure orchestrator maps caller name -> method, translates options,
% runs bst_eigen per input, and rejects unknown callers without throwing.
% SKIPs if the live protocol lacks an LBO eigen node + a constrained results file on its surface.
    nFail = 0; chk = @i_chk;
    [ef, surf, nv] = i_find_eigen('Laplace-Beltrami');
    if isempty(ef); fprintf('SKIP test_process_eigen: no LBO eigen node.\n'); return; end
    rf = i_find_results(surf, nv);
    if isempty(rf); fprintf('SKIP test_process_eigen: no constrained results on surface.\n'); return; end

    % Synthetic sProcess mimicking process_eigenspectrum (str2func name need not exist on disk)
    sProcess = struct();
    sProcess.Function = str2func('process_eigenspectrum');
    sProcess.options.variant.Value     = 'Laplace-Beltrami';
    sProcess.options.measure.Value     = 'power';
    sProcess.options.timewindow.Value  = [];
    sProcess.options.win_length.Value  = {[], 's', []};
    sProcess.options.win_overlap.Value = {50, '%', 1};
    sProcess.options.win_std.Value     = 0;
    sInputs = bst_process('GetInputStruct', rf);

    Out = process_eigen('Run', sProcess, sInputs);
    ok = iscell(Out) && ~isempty(Out) && ischar(Out{1});
    nFail = nFail + chk('process_eigen returns a saved node', ok);
    if ok
        R = in_bst_timefreq(Out{1}, 'Method', 'Freqs');
        nFail = nFail + chk('node Method == eigenspectrum', strcmp(R.Method, 'eigenspectrum'));
        nFail = nFail + chk('node has Freqs', ~isempty(R.Freqs));
        try, file_delete(file_fullpath(Out{1}), 1); db_reload_studies(sInputs.iStudy); catch; end
    end

    % Unknown caller -> empty output, no throw
    sBad = sProcess; sBad.Function = str2func('process_eigenbogus');
    ob = process_eigen('Run', sBad, sInputs);
    nFail = nFail + chk('unknown caller -> empty output', isempty(ob));

    fprintf('\n==== test_process_eigen: %d failed ====\n', nFail);
    if nFail > 0; error('test_process_eigen FAILED'); end
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
