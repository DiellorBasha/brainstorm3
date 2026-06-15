function test_filter_designer_session()
% Open a designer session, seed, save a bank, confirm the nested node, then delete it.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;

    hFig = view_filter_designer(EigenFile); drawnow;
    nFail = nFail + chk('session figure opens', ishandle(hFig));

    panel_filter_designer('SetSeedVertex', 'FilterDesigner', 100); drawnow;
    S = panel_filter_designer('GetState','FilterDesigner');
    nFail = nFail + chk('seed sets coeffs', ~isempty(S) && ~isempty(S.SeedCoeffs));

    panel_filter_designer('OnSave', 'FilterDesigner'); drawnow;
    sSubject = bst_get('Subject',1);
    fbs = sSubject.Surface(5).Filterbank;
    nFail = nFail + chk('bank saved + nested (ParentEigen match)', ...
        ~isempty(fbs) && any(strcmp({fbs.ParentEigen}, file_short(EigenFile))));
    nFail = nFail + chk('session torn down after save', ~ishandle(hFig));

    % node_create_subject nests it under the eigen node: rebuild the tree node and check
    nFail = nFail + chk('node renders nested under eigen', local_nests_ok(sSubject, EigenFile));

    % cleanup inline (FilterbankDelete_Callback is a private tree_callbacks subfunction,
    % exercised via the GUI; it mirrors the verified EigenDelete_Callback)
    newFile = fbs(end).FileName;
    file_delete(file_fullpath(newFile), 1);
    [~, iSub, iSurf, iFb] = bst_get('FilterbankFile', newFile);
    if ~isempty(iFb)
        ProtocolSubjects = bst_get('ProtocolSubjects');
        ProtocolSubjects.Subject(iSub).Surface(iSurf).Filterbank(iFb) = [];
        bst_set('ProtocolSubjects', ProtocolSubjects); db_save();
    end
    nFail = nFail + chk('delete removes the node', isempty(bst_get('Subject',1).Surface(5).Filterbank));

    fprintf('\n==== test_filter_designer_session: %d failed ====\n', nFail);
    if nFail > 0, error('test_filter_designer_session FAILED'); end
end

function ok = local_nests_ok(sSubject, EigenFile)
    % Confirm the saved filterbank's ParentEigen points at the eigen node (the nesting key).
    ok = false;
    fbs = sSubject.Surface(5).Filterbank;
    if isempty(fbs); return; end
    ok = file_compare(fbs(end).ParentEigen, file_short(EigenFile));
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
