function preventad_import_batch(BidsDir, LogFile)
% PREVENTAD_IMPORT_BATCH: Import all FreeSurfer PREVENT-AD subjects into 'preventad'.
%
% Auto-discovers every subject that has FreeSurfer anatomy (icosphere/ico5 is
% FreeSurfer-only), then runs preventad_import on each one that is NOT already in
% the protocol. Designed for an unattended overnight run:
%   - RESUMABLE  : subjects already present in the protocol are skipped, so the
%                  batch can be stopped and restarted without redoing work.
%   - ROBUST     : each subject runs in try/catch, so one failure does not abort
%                  the batch; the failure is logged and the run continues.
%   - MONITORABLE: timestamped per-subject progress is appended to LogFile
%                  (readable from outside MATLAB while the batch runs).
%
% USAGE:  preventad_import_batch()                 % default BidsDir + default log
%         preventad_import_batch(BidsDir)
%         preventad_import_batch(BidsDir, LogFile)
%
% INPUTS:
%     BidsDir : PREVENT-AD BIDS MEG root (default: the SpikeData-2 dataset path).
%     LogFile : Progress log path (default: dev/preventad_import_batch.log next to
%               this script).
%
% Author: Diellor Basha, 2026

    if (nargin < 1) || isempty(BidsDir)
        BidsDir = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/meg';
    end
    if (nargin < 2) || isempty(LogFile)
        LogFile = bst_fullfile(bst_fileparts(mfilename('fullpath')), 'preventad_import_batch.log');
    end
    if ~file_exist(BidsDir)
        error('BIDS directory not found: %s', BidsDir);
    end

    % --- discover FreeSurfer subjects (icosphere/ico5 is FreeSurfer-only) ---
    fsDir   = bst_fullfile(BidsDir, 'derivatives', 'freesurfer');
    fsSubj  = dir(bst_fullfile(fsDir, 'sub-*'));
    fsSubj  = fsSubj([fsSubj.isdir]);
    SubjectNames = sort({fsSubj.name});
    nTotal = numel(SubjectNames);

    % --- names already present in the protocol (skip these) ---
    PS = bst_get('ProtocolSubjects');
    existing = {};
    if ~isempty(PS) && ~isempty(PS.Subject)
        existing = {PS.Subject.Name};
    end

    log_line(LogFile, sprintf('==== BATCH START: %d FreeSurfer subjects, %d already in protocol ====', ...
        nTotal, numel(intersect(SubjectNames, existing))));

    nDone = 0; nSkip = 0; nFail = 0;
    tBatch = tic;
    for i = 1:nTotal
        name = SubjectNames{i};
        if ismember(name, existing)
            nSkip = nSkip + 1;
            log_line(LogFile, sprintf('[%d/%d] SKIP (already imported): %s', i, nTotal, name));
            continue;
        end
        log_line(LogFile, sprintf('[%d/%d] START: %s', i, nTotal, name));
        tSubj = tic;
        try
            preventad_import(BidsDir, name, false);   % DoSnapshots=false: headless-safe (no hanging 3D renders)
            nDone = nDone + 1;
            log_line(LogFile, sprintf('[%d/%d] DONE: %s (%.1f min)', i, nTotal, name, toc(tSubj)/60));
        catch ME
            nFail = nFail + 1;
            log_line(LogFile, sprintf('[%d/%d] FAIL: %s (%.1f min) - %s', i, nTotal, name, toc(tSubj)/60, ME.message));
            for k = 1:numel(ME.stack)
                log_line(LogFile, sprintf('        at %s (line %d)', ME.stack(k).name, ME.stack(k).line));
            end
        end
    end

    log_line(LogFile, sprintf('==== BATCH END: %d imported, %d skipped, %d failed, total %.1f h ====', ...
        nDone, nSkip, nFail, toc(tBatch)/3600));
end


%% ===== append one timestamped line to the log (and echo to command window) =====
function log_line(LogFile, msg)
    stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    line  = sprintf('%s  %s', stamp, msg);
    disp(line);
    fid = fopen(LogFile, 'a');
    if fid > 0
        fprintf(fid, '%s\n', line);
        fclose(fid);
    end
end
