function test_db_add_wavelet()
% Round-trip: save a wavelet under an eigen node, resolve it, reload it, delete it.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    sSubject = bst_get('Subject', 1);
    iSurf = 5;
    assert(isfield(sSubject.Surface(iSurf),'Eigen') && ~isempty(sSubject.Surface(iSurf).Eigen), ...
        'No eigen node on Subject01 surface 5; compute one first.');
    EigenFile = sSubject.Surface(iSurf).Eigen(1).FileName;

    wavelet = struct('Kernel','mexhat','Params',struct('t',0.01),'Direction',[1 0 0], ...
                     'Chirality',0,'Axis',[0 0 1]);
    opts    = struct('N',4,'Spacing','geometric','LambdaRange',[1 256],'Chiralities',[]);
    w = db_template('waveletmat');
    w.ParentEigen = file_short(EigenFile);
    w.Variant     = 'Dirac';
    w.Tiles       = bst_filterbank_tiles(wavelet, opts);
    w.Tiling      = struct('Wavelet', wavelet, 'Opts', opts);

    iW = db_add_wavelet(1, EigenFile, w, 'TEST wavelet');
    nFail = nFail + chk('returns an index', ~isempty(iW));

    sSubject = bst_get('Subject', 1);
    ws = sSubject.Surface(iSurf).Wavelet;
    nFail = nFail + chk('appears in Surface.Wavelet', ~isempty(ws));
    nFail = nFail + chk('ParentEigen matches', any(strcmp({ws.ParentEigen}, file_short(EigenFile))));

    newFile = ws(end).FileName;
    [~,~,iSurf2,iW2] = bst_get('WaveletFile', newFile);
    nFail = nFail + chk('accessor resolves it', ~isempty(iW2) && iSurf2==iSurf);
    R = load(file_fullpath(newFile));
    nFail = nFail + chk('reloaded Tiles count', numel(R.Tiles)==numel(w.Tiles));
    nFail = nFail + chk('reloaded Tiling.Wavelet', strcmp(R.Tiling.Wavelet.Kernel,'mexhat'));

    file_delete(file_fullpath(newFile), 1);
    sSubject.Surface(iSurf).Wavelet(end) = [];
    ProtocolSubjects = bst_get('ProtocolSubjects');
    ProtocolSubjects.Subject(1) = sSubject;
    bst_set('ProtocolSubjects', ProtocolSubjects);
    db_save();

    fprintf('\n==== test_db_add_wavelet: %d failed ====\n', nFail);
    if nFail > 0, error('test_db_add_wavelet FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
