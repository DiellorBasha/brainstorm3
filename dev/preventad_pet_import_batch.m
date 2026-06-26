function preventad_pet_import_batch(BidsPetDir, LogFile, Opts)
% PREVENTAD_PET_IMPORT_BATCH: Import + process PREVENT-AD PET for every subject that
% ALREADY EXISTS in the 'preventad' protocol (i.e. has MEG/FreeSurfer anatomy),
% registering each subject's PET tracers to that existing anatomy.
%
% Auto-discovers PET subjects (BIDS sub-* under BidsPetDir), and for each one that
% is present in the protocol, runs preventad_pet_import. Built for an unattended run:
%   - MATCHED ONLY : PET subjects with no existing anatomy are skipped (logged) -
%                    PET registers to the MEG anatomy, so unmatched subjects are out.
%   - RESUMABLE    : preventad_pet_import skips tracers whose SUVR already exists,
%                    so the batch can be stopped/restarted without redoing work.
%   - ROBUST       : each subject runs in try/catch; a failure is logged, run continues.
%   - MONITORABLE  : timestamped per-subject progress appended to LogFile.
%
% USAGE:  preventad_pet_import_batch()
%         preventad_pet_import_batch(BidsPetDir)
%         preventad_pet_import_batch(BidsPetDir, LogFile)
%         preventad_pet_import_batch(BidsPetDir, LogFile, Opts)   % Opts -> preventad_pet_import
%
% INPUTS:
%   BidsPetDir : PREVENT-AD PET BIDS root (default: SpikeData-2 pet dataset).
%   LogFile    : Progress log path (default: dev/preventad_pet_import_batch.log).
%   Opts       : Options struct forwarded to preventad_pet_import (atlas/ROI/etc.).
%
% Author: Diellor Basha, 2026

    if (nargin < 1) || isempty(BidsPetDir)
        BidsPetDir = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet';
    end
    if (nargin < 2) || isempty(LogFile)
        LogFile = bst_fullfile(bst_fileparts(mfilename('fullpath')), 'preventad_pet_import_batch.log');
    end
    if (nargin < 3) || isempty(Opts), Opts = struct(); end
    if ~file_exist(BidsPetDir)
        error('PET BIDS directory not found: %s', BidsPetDir);
    end

    % ---- select existing protocol ----
    ProtocolName = 'preventad';
    if ~brainstorm('status'), brainstorm nogui; end
    iProtocol = bst_get('Protocol', ProtocolName);
    if isempty(iProtocol)
        error('Unknown protocol: %s (run the MEG import first).', ProtocolName);
    end
    if (iProtocol ~= bst_get('iProtocol'))
        gui_brainstorm('SetCurrentProtocol', iProtocol);
    end

    % ---- discover PET subjects (skip AppleDouble '._*') ----
    d = dir(bst_fullfile(BidsPetDir, 'sub-*'));
    d = d([d.isdir]);
    names = {d.name};
    names = names(~strncmp(names, '._', 2));
    SubjectNames = sort(names);
    nTotal = numel(SubjectNames);

    log_line(LogFile, sprintf('==== PET BATCH START: %d PET subjects discovered ====', nTotal));

    nDone = 0; nSkip = 0; nFail = 0;
    tBatch = tic;
    for i = 1:nTotal
        name = SubjectNames{i};
        % match to an EXISTING subject (PET registers to the MEG anatomy)
        [sSubject, iSubject] = bst_get('Subject', name);
        if isempty(iSubject) || isempty(sSubject) || isempty(sSubject.Anatomy)
            nSkip = nSkip + 1;
            log_line(LogFile, sprintf('[%d/%d] SKIP (no matching anatomy): %s', i, nTotal, name));
            continue;
        end
        log_line(LogFile, sprintf('[%d/%d] START: %s', i, nTotal, name));
        tSubj = tic;
        try
            R = preventad_pet_import(BidsPetDir, name, Opts);
            nTrac = sum(~[R.Skipped]); nSk = sum([R.Skipped]);
            nDone = nDone + 1;
            log_line(LogFile, sprintf('[%d/%d] DONE: %s (%d tracers, %d already-done, %.1f min)', ...
                i, nTotal, name, nTrac, nSk, toc(tSubj)/60));
        catch ME
            nFail = nFail + 1;
            log_line(LogFile, sprintf('[%d/%d] FAIL: %s (%.1f min) - %s', i, nTotal, name, toc(tSubj)/60, ME.message));
            for k = 1:numel(ME.stack)
                log_line(LogFile, sprintf('        at %s (line %d)', ME.stack(k).name, ME.stack(k).line));
            end
        end
    end

    log_line(LogFile, sprintf('==== PET BATCH END: %d subjects done, %d skipped(no anat), %d failed, total %.1f h ====', ...
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
