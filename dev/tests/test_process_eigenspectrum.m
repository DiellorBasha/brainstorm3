function test_process_eigenspectrum()
% TEST_PROCESS_EIGENSPECTRUM: the registered process advertises the right contract and runs
% end-to-end through process_eigen -> bst_eigen, producing a timefreq_eigenspectrum node.
% End-to-end portion SKIPs if the live protocol lacks the fixture.
    nFail = 0; chk = @i_chk;

    % ---- GetDescription contract ----
    sProcess = process_eigenspectrum('GetDescription');
    nFail = nFail + chk('Comment',     strcmp(sProcess.Comment, 'Eigenspectrum (spatial PSD)'));
    nFail = nFail + chk('SubGroup',    strcmp(sProcess.SubGroup, 'Frequency'));
    nFail = nFail + chk('Index 484',   isequal(sProcess.Index, 484));
    nFail = nFail + chk('InputTypes',  isequal(sProcess.InputTypes, {'results'}));
    nFail = nFail + chk('OutputTypes', isequal(sProcess.OutputTypes, {'timefreq'}));
    nFail = nFail + chk('variant default LBO', strcmp(sProcess.options.variant.Value, 'Laplace-Beltrami'));
    nFail = nFail + chk('measure default power', strcmp(sProcess.options.measure.Value, 'power'));
    nFail = nFail + chk('win_length present', isfield(sProcess.options, 'win_length'));
    nFail = nFail + chk('win_overlap present', isfield(sProcess.options, 'win_overlap'));
    nFail = nFail + chk('win_std present', isfield(sProcess.options, 'win_std'));

    % ---- end-to-end via the registered process ----
    [ef, surf, nv] = i_find_eigen('Laplace-Beltrami');
    if isempty(ef)
        fprintf('SKIP e2e: no LBO eigen node.\n');
    else
        rf = i_find_results(surf, nv);
        if isempty(rf)
            fprintf('SKIP e2e: no constrained results on surface.\n');
        else
            sProcess.Function = @process_eigenspectrum;
            sInputs = bst_process('GetInputStruct', rf);
            Out = process_eigenspectrum('Run', sProcess, sInputs);
            ok = iscell(Out) && ~isempty(Out) && ischar(Out{1});
            nFail = nFail + chk('Run produces a node', ok);
            if ok
                nFail = nFail + chk('node prefix timefreq_eigenspectrum', ...
                    ~isempty(strfind(Out{1}, 'timefreq_eigenspectrum'))); %#ok<STREMP>
                R = in_bst_timefreq(Out{1}, 0, 'Method');
                nFail = nFail + chk('node Method == eigenspectrum', strcmp(R.Method, 'eigenspectrum'));
                try, file_delete(file_fullpath(Out{1}), 1); db_reload_studies(sInputs.iStudy); catch; end
            end
        end
    end

    fprintf('\n==== test_process_eigenspectrum: %d failed ====\n', nFail);
    if nFail > 0; error('test_process_eigenspectrum FAILED'); end
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
