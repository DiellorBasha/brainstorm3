function test_wavelet_designer_session()
% Open a designer session, seed, save a wavelet, confirm the nested node, then delete it.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;
    nWStart = numel(bst_get('Subject',1).Surface(5).Wavelet);   % isolation: count before

    hFig = view_wavelet_designer(EigenFile); drawnow;
    nFail = nFail + chk('session figure opens', ishandle(hFig));

    panel_wavelet_designer('SetSeedVertex', 'WaveletDesigner', 100); drawnow;
    S = panel_wavelet_designer('GetState','WaveletDesigner');
    nFail = nFail + chk('seed sets coeffs', ~isempty(S) && ~isempty(S.SeedCoeffs));

    panel_wavelet_designer('OnSave', 'WaveletDesigner'); drawnow;
    sSubject = bst_get('Subject',1);
    ws = sSubject.Surface(5).Wavelet;
    nFail = nFail + chk('wavelet saved + nested (ParentEigen match)', ...
        ~isempty(ws) && any(strcmp({ws.ParentEigen}, file_short(EigenFile))));
    nFail = nFail + chk('session torn down after save', ~ishandle(hFig));

    % node_create_subject nests it under the eigen node: confirm the nesting key
    nFail = nFail + chk('node renders nested under eigen', local_nests_ok(sSubject, EigenFile));

    % cleanup inline (WaveletDelete_Callback is a private tree_callbacks subfunction,
    % exercised via the GUI; it mirrors the verified EigenDelete_Callback)
    newFile = ws(end).FileName;
    file_delete(file_fullpath(newFile), 1);
    [~, iSub, iSurf, iW] = bst_get('WaveletFile', newFile);
    if ~isempty(iW)
        ProtocolSubjects = bst_get('ProtocolSubjects');
        ProtocolSubjects.Subject(iSub).Surface(iSurf).Wavelet(iW) = [];
        bst_set('ProtocolSubjects', ProtocolSubjects); db_save();
    end
    nWEnd = numel(bst_get('Subject',1).Surface(5).Wavelet);
    nFail = nFail + chk('delete removes the test node (count restored)', nWEnd == nWStart);
    nFail = nFail + chk('test node file gone', ~file_exist(file_fullpath(newFile)));

    fprintf('\n==== test_wavelet_designer_session: %d failed ====\n', nFail);
    if nFail > 0, error('test_wavelet_designer_session FAILED'); end
end

function ok = local_nests_ok(sSubject, EigenFile)
    % Confirm the saved wavelet's ParentEigen points at the eigen node (the nesting key).
    ok = false;
    ws = sSubject.Surface(5).Wavelet;
    if isempty(ws); return; end
    ok = file_compare(ws(end).ParentEigen, file_short(EigenFile));
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
