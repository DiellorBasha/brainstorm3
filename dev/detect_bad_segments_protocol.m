function report = detect_bad_segments_protocol(ProtocolZip, OutDir, BstDbDir)
% DETECT_BAD_SEGMENTS_PROTOCOL  Brainstorm's bad-segment detector on the resting MEG of one exported
%                               per-subject protocol; each recording's events saved for the converter.
%
%   report = detect_bad_segments_protocol(ProtocolZip, OutDir, BstDbDir)
%
% The same step as the multimodal worker's NSP_BAD_SEGMENTS (preventad_multimodal_subject.m), standalone
% so it can run on a protocol that is already built (locally, headless):
%   process_evt_detect_badsegment  on every task-rest raw recording (MEG; 1–7 Hz eye/movement and
%                                  40–240 Hz muscle; sensitivity 3 — preventad_import.m's settings)
%   process_evt_rename             1-7Hz, 40-240Hz → bad_1-7Hz, bad_40-240Hz (Brainstorm's "bad" rule)
% A recording that already carries bad* events is left alone.
%
% OUTPUT (OutDir):
%   <condition>_events.mat   the recording's full F.events after detection (the struct the converter
%                            reads when it writes timeseries/<rec>_events), plus 'condition' and 'file'
%   bad_segments.json        per recording: segments, seconds, whether they were detected or kept
%
% Run headless with the Brainstorm home isolated (never the user's own ~/.brainstorm):
%   JAVA_TOOL_OPTIONS=-Duser.home=<dir> matlab -batch "addpath('<bst>'); addpath('<bst>/dev'); ..."
%
% Authors: Diellor Basha, 2026

    if exist(OutDir, 'dir') ~= 7, mkdir(OutDir); end
    if exist(BstDbDir, 'dir') ~= 7, mkdir(BstDbDir); end
    brainstorm setpath;
    bst_user_dir = fullfile(char(java.lang.System.getProperty('user.home')), '.brainstorm');
    if exist(bst_user_dir, 'dir') ~= 7, mkdir(bst_user_dir); end
    iProtocol             = 0; %#ok<NASGU>
    ProtocolsListInfo     = repmat(db_template('ProtocolInfo'), 0);     %#ok<NASGU>
    ProtocolsListSubjects = repmat(db_template('ProtocolSubjects'), 0); %#ok<NASGU>
    ProtocolsListStudies  = repmat(db_template('ProtocolStudies'), 0);  %#ok<NASGU>
    BrainStormDbDir       = BstDbDir; %#ok<NASGU>
    DbVersion             = 5.03; %#ok<NASGU>
    save(fullfile(bst_user_dir, 'brainstorm.mat'), 'iProtocol', 'ProtocolsListInfo', ...
         'ProtocolsListSubjects', 'ProtocolsListStudies', 'BrainStormDbDir', 'DbVersion');
    if ~brainstorm('status'); brainstorm server; end
    report = struct('condition', {}, 'segments', {}, 'seconds', {}, 'detected', {});
    try
        import_protocol(ProtocolZip);
        sProt = bst_get('ProtocolSubjects');
        for iSubject = 1:numel(sProt.Subject)
            sSubject = bst_get('Subject', iSubject);
            sStudies = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject');
            for i = 1:numel(sStudies)
                for k = 1:numel(sStudies(i).Data)
                    f = sStudies(i).Data(k).FileName;
                    if isempty(strfind(f, 'task-rest')) || ~strcmpi(sStudies(i).Data(k).DataType, 'raw'), continue; end
                    ev = local_events(f);
                    has = any(strncmpi({ev.label}, 'bad', 3));
                    if ~has
                        bst_process('CallProcess', 'process_evt_detect_badsegment', {f}, [], ...
                            'timewindow', [], 'sensortypes', 'MEG', 'threshold', 3, 'isLowFreq', 1, 'isHighFreq', 1);
                        bst_process('CallProcess', 'process_evt_rename', {f}, [], ...
                            'src', '1-7Hz, 40-240Hz', 'dest', 'bad_1-7Hz, bad_40-240Hz');
                        ev = local_events(f);
                    end
                    bad = ev(strncmpi({ev.label}, 'bad', 3));
                    n = 0; sec = 0;
                    for b = 1:numel(bad)
                        n = n + size(bad(b).times, 2);
                        if size(bad(b).times, 1) == 2, sec = sec + sum(diff(bad(b).times, 1, 1)); end
                    end
                    condition = sStudies(i).Condition{1};
                    if strncmp(condition, '@raw', 4), condition = condition(5:end); end
                    events = ev; file = f; %#ok<NASGU>
                    save(fullfile(OutDir, [condition '_events.mat']), 'events', 'condition', 'file', '-v7');
                    report(end+1) = struct('condition', condition, 'segments', n, 'seconds', sec, 'detected', ~has); %#ok<AGROW>
                    fprintf('bad segments: %s -> %d segment(s), %.1f s\n', condition, n, sec);
                end
            end
        end
        fid = fopen(fullfile(OutDir, 'bad_segments.json'), 'w');
        fprintf(fid, '%s', jsonencode(report));
        fclose(fid);
        brainstorm stop;
    catch ME
        try brainstorm stop; catch; end
        rethrow(ME);
    end
end

function ev = local_events(dataFile)
    D = in_bst_data(dataFile, 'F');
    if isfield(D.F, 'events') && ~isempty(D.F.events)
        ev = D.F.events;
    else
        ev = struct('label', {}, 'times', {});
    end
end
