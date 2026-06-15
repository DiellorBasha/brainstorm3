function test_db_add_filterbank()
% Round-trip: save a filterbank under an eigen node, resolve it, reload it, delete it.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    sSubject = bst_get('Subject', 1);
    iSurf = 5;
    assert(isfield(sSubject.Surface(iSurf),'Eigen') && ~isempty(sSubject.Surface(iSurf).Eigen), ...
        'No eigen node on Subject01 surface 5; compute one first.');
    EigenFile = sSubject.Surface(iSurf).Eigen(1).FileName;

    base = struct('Kernel','mexhat','Params',struct(),'Direction',[1 0 0], ...
                  'Chirality',0,'Axis',[0 0 1],'N',4,'Spacing','geometric', ...
                  'LambdaRange',[1 256],'Chiralities',0);
    fb = db_template('filterbankmat');
    fb.ParentEigen = file_short(EigenFile);
    fb.Variant     = 'Dirac';
    fb.Tiles       = bst_filterbank_tiles(base);
    fb.Tiling      = base;

    iFb = db_add_filterbank(1, EigenFile, fb, 'TEST filterbank');
    nFail = nFail + chk('returns an index', ~isempty(iFb));

    % registered under the same surface, keyed by ParentEigen
    sSubject = bst_get('Subject', 1);
    fbs = sSubject.Surface(iSurf).Filterbank;
    nFail = nFail + chk('appears in Surface.Filterbank', ~isempty(fbs));
    nFail = nFail + chk('ParentEigen matches', any(strcmp({fbs.ParentEigen}, file_short(EigenFile))));

    % resolvable via accessor + reloads to an identical recipe bank
    newFile = fbs(end).FileName;
    [s2,iSub2,iSurf2,iFb2] = bst_get('FilterbankFile', newFile);
    nFail = nFail + chk('accessor resolves it', ~isempty(iFb2) && iSurf2==iSurf);
    R = load(file_fullpath(newFile));
    nFail = nFail + chk('reloaded Tiles count', numel(R.Tiles)==numel(fb.Tiles));
    nFail = nFail + chk('reloaded Variant', strcmp(R.Variant,'Dirac'));

    % cleanup: delete the test node + file (proper DB path)
    file_delete(file_fullpath(newFile), 1);
    sSubject.Surface(iSurf).Filterbank(end) = [];
    ProtocolSubjects = bst_get('ProtocolSubjects');
    ProtocolSubjects.Subject(1) = sSubject;
    bst_set('ProtocolSubjects', ProtocolSubjects);
    db_save();

    fprintf('\n==== test_db_add_filterbank: %d failed ====\n', nFail);
    if nFail > 0, error('test_db_add_filterbank FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
