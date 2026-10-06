function preventad_subject(BidsDir, OutputDir, SubjectLabel, Module, varargin)
% PREVENTAD_SUBJECT  Per-subject PREVENT-AD MEG pipeline for nsp/HPC/container use.
%
% A single-subject, self-contained, throwaway-protocol adaptation of
% dev/preventad_import.m, built on the same container contract as
% brainstorm-container/scripts/bst_single_subject.m so it drops straight into the
% nsp per-subject pattern (one SLURM array task = one subject).
%
% Differences from bst_single_subject (the PREVENT-AD customizations):
%   1. Anatomy is downsampled with the FreeSurfer ICOSPHERE method at ICO5
%      (20484 vtx) — HARD REQUIREMENT for the differential/spectral engine
%      (reducepatch breaks the sphere-registered correspondence). Uses
%      downsamplemethod='icosphere' + icolevel=ICO_LEVEL, NOT nvertices.
%   2. Source stage computes the nxr manifold backbone (process_tess_manifold)
%      and BOTH inverses on one overlapping-spheres head model:
%        - standard dSPM        (process_inverse_2018, unconstrained)
%        - custom Dirac dSPM     (process_inverse_dirac)
%      => requires the nxr-compute plugin (and SPM12 for anatomy).
%   3. Noise covariance comes from the within-subject task-noise recording.
%
% USAGE (matches bst_single_subject so the container can call it the same way):
%   preventad_subject(BidsDir, OutputDir, SubjectLabel, Module)
%   preventad_subject(BidsDir, OutputDir, SubjectLabel, Module, 'Key', Value, ...)
%
% REQUIRED:
%   BidsDir      - BIDS MEG root (the per-subject clone; must contain sub-<label>/
%                  and derivatives/freesurfer/sub-<label>/)
%   OutputDir    - Where to write sub-<label>_brainstorm.zip
%   SubjectLabel - Subject label WITHOUT 'sub-' prefix (e.g., 'MTL0002')
%   Module       - Stop position: 'import' | 'preprocess' | 'source' | 'timefreq'
%
% OPTIONAL key-value:
%   'BstDir'    - brainstorm3 source tree (the DEV FORK; default: auto-detect)
%   'BstDbDir'  - throwaway protocol DB (default: $SLURM_TMPDIR/brainstorm_db)
%   'IcoLevel'  - icosphere level: 'ico3'|'ico4'|'ico5'|'ico6' (default 'ico5')
%   'NVertices' - accepted for bst_single_subject API parity; IGNORED on the
%                 icosphere path (kept so the container wrapper can pass it).
%
% OUTPUT:  <OutputDir>/sub-<SubjectLabel>_brainstorm.zip
%
% CONTAINER INTEGRATION (one-line change): in the brainstorm-container entrypoint,
% point participant mode at this function instead of bst_single_subject, e.g.
%   addpath(BstDir); addpath(fullfile(BstDir,'dev'));
%   preventad_subject(bids_dir, out_dir, label, module, 'BstDir',BstDir, 'BstDbDir',db);
% (preventad_subject lives in the fork at dev/, so binding --bst-dir <fork> exposes it.)
%
% Author: Diellor Basha, 2026  (scaffold: bst_single_subject; chain: preventad_import)

%% ===== PARSE INPUTS =====
p = inputParser;
addRequired(p, 'BidsDir', @ischar);
addRequired(p, 'OutputDir', @ischar);
addRequired(p, 'SubjectLabel', @ischar);
addRequired(p, 'Module', @(x) ismember(x, {'import','preprocess','source','timefreq'}));
addParameter(p, 'BstDir', '', @ischar);
addParameter(p, 'BstDbDir', '', @ischar);
addParameter(p, 'IcoLevel', 'ico5', @ischar);
if isdeployed
    addParameter(p, 'NVertices', 15000, @(x) isnumeric(x) || ischar(x));
else
    addParameter(p, 'NVertices', 15000, @isnumeric);
end
parse(p, BidsDir, OutputDir, SubjectLabel, Module, varargin{:});
opts = p.Results;

SubjectName  = ['sub-' opts.SubjectLabel];
ProtocolName = ['nsp_' SubjectName];

fprintf('\n================================================================\n');
fprintf(' PREVENTAD_SUBJECT — Per-subject PREVENT-AD MEG pipeline\n');
fprintf('================================================================\n');
fprintf(' Subject:    %s\n', SubjectName);
fprintf(' Module:     %s (stop position)\n', opts.Module);
fprintf(' IcoLevel:   %s (icosphere; nvertices ignored)\n', opts.IcoLevel);
fprintf(' BIDS dir:   %s\n', opts.BidsDir);
fprintf(' Output dir: %s\n', opts.OutputDir);
fprintf('================================================================\n\n');

%% ===== RESOLVE PATHS (expand ~ — Java IO does not handle tilde) =====
opts.BidsDir   = expand_tilde(opts.BidsDir);
opts.OutputDir = expand_tilde(opts.OutputDir);
opts.BstDir    = expand_tilde(opts.BstDir);
opts.BstDbDir  = expand_tilde(opts.BstDbDir);

% --- Brainstorm source tree (the DEV FORK) ---
BstDir = opts.BstDir;
if isempty(BstDir)
    if exist('brainstorm', 'file') == 2
        BstDir = fileparts(which('brainstorm'));
    else
        error('BstDir not specified and brainstorm3 not on MATLAB path');
    end
elseif ~isdeployed
    addpath(BstDir);
    if exist(fullfile(BstDir,'dev'), 'dir'); addpath(fullfile(BstDir,'dev')); end
end

% --- throwaway DB dir (node-local) ---
BstDbDir = opts.BstDbDir;
if isempty(BstDbDir)
    slurm_tmpdir = getenv('SLURM_TMPDIR');
    if ~isempty(slurm_tmpdir)
        BstDbDir = fullfile(slurm_tmpdir, 'brainstorm_db');
    else
        BstDbDir = fullfile(tempdir, 'brainstorm_db');
    end
end
if ~exist(BstDbDir, 'dir');   mkdir(BstDbDir);   end
if ~exist(opts.OutputDir,'dir'); mkdir(opts.OutputDir); end

%% ===== LOGGING =====
LogFile = fullfile(opts.OutputDir, [SubjectName '_log.txt']);
diary(LogFile);
fprintf('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% ===== MODULES TO RUN (module = stop position) =====
MODULE_ORDER  = {'import', 'preprocess', 'source', 'timefreq'};
target_idx    = find(strcmp(MODULE_ORDER, opts.Module));
modules_to_run = MODULE_ORDER(1:target_idx);
fprintf('Pipeline: %s\n', strjoin(modules_to_run, ' -> '));

%% ===== INIT BRAINSTORM (headless server) + throwaway protocol config =====
brainstorm setpath;
bst_user_dir = fullfile(char(java.lang.System.getProperty('user.home')), '.brainstorm');
if ~exist(bst_user_dir, 'dir'); mkdir(bst_user_dir); end
bst_cfg_file          = fullfile(bst_user_dir, 'brainstorm.mat'); %#ok<NASGU>
iProtocol             = 0; %#ok<NASGU>
ProtocolsListInfo     = repmat(db_template('ProtocolInfo'), 0); %#ok<NASGU>
ProtocolsListSubjects = repmat(db_template('ProtocolSubjects'), 0); %#ok<NASGU>
ProtocolsListStudies  = repmat(db_template('ProtocolStudies'), 0); %#ok<NASGU>
BrainStormDbDir       = BstDbDir; %#ok<NASGU>
DbVersion             = 5.03; %#ok<NASGU>
save(fullfile(bst_user_dir,'brainstorm.mat'), 'iProtocol', 'ProtocolsListInfo', ...
     'ProtocolsListSubjects', 'ProtocolsListStudies', 'BrainStormDbDir', 'DbVersion');
if ~brainstorm('status'); brainstorm server; end
fprintf('Brainstorm server started. DB = %s\n', BstDbDir);

%% ===== CREATE THROWAWAY PROTOCOL =====
iExisting = bst_get('Protocol', ProtocolName);
if ~isempty(iExisting); gui_brainstorm('DeleteProtocol', ProtocolName); end
protocolDir = fullfile(BstDbDir, ProtocolName);
if exist(protocolDir, 'dir'); rmdir(protocolDir, 's'); end
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);   % UseDefaultAnat=0, UseDefaultChannel=0
bst_report('Start');
fprintf('Protocol created: %s\n', ProtocolName);

%% ===== RUN MODULES =====
sFilesRaw = []; sFilesBand = []; sFilesRest = []; sSrc = []; %#ok<NASGU>
try
    for iMod = 1:numel(modules_to_run)
        mod_name = modules_to_run{iMod};
        fprintf('\n======== Module: %s (%d/%d) ========\n', upper(mod_name), iMod, numel(modules_to_run));
        switch mod_name
            case 'import'
                sFilesRaw  = pa_import(opts.BidsDir, SubjectName, opts.IcoLevel);
            case 'preprocess'
                [sFilesRest, sFilesBand] = pa_preprocess(sFilesRaw);
            case 'source'
                sSrc       = pa_source(sFilesRest, sFilesBand); %#ok<NASGU>
            case 'timefreq'
                pa_power(sFilesRest);
        end
    end
    fprintf('\n--- All modules completed ---\n');
catch ME
    fprintf('\nERROR in module "%s": %s\n', mod_name, ME.message);
    for k = 1:numel(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    R = bst_report('Save');
    if ~isempty(R); try, copyfile(R, fullfile(opts.OutputDir,[SubjectName '_report_error.html'])); end; end %#ok<TRYNC>
    brainstorm stop; diary off; rethrow(ME);
end

%% ===== EXPORT SUBJECT AS SELF-CONTAINED .zip =====
[~, iSubject] = bst_get('Subject', SubjectName);
if isempty(iSubject) || iSubject == 0
    error('Could not find subject %s in protocol for export', SubjectName);
end
iProtocol  = bst_get('iProtocol');
ExportZip  = fullfile(opts.OutputDir, [SubjectName '_brainstorm.zip']);
export_protocol(iProtocol, iSubject, ExportZip);
zi = dir(ExportZip);
if isempty(zi) || zi.bytes == 0
    error('Export failed — zip missing/empty: %s', ExportZip);
end
fprintf('Export complete: %s (%.1f MB)\n', ExportZip, zi.bytes/1048576);

%% ===== SAVE REPORT + CLEANUP =====
R = bst_report('Save');
if ~isempty(R); try, copyfile(R, fullfile(opts.OutputDir,[SubjectName '_report.html'])); end; end %#ok<TRYNC>
brainstorm stop;
fprintf('Finished: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
diary off;
end


%% ########################################################################
%% MODULE: IMPORT  (icosphere/ico5 FreeSurfer anatomy)
%% ########################################################################
function sFilesRaw = pa_import(BidsDir, SubjectName, IcoLevel)
% selectsubj restricts to the single subject; icosphere/ico5 is the hard
% requirement; anatregister left at default so SPM12 handles MRI/MNI as in
% preventad_import (SPM12 must be installed as a Brainstorm plugin).
sFilesRaw = bst_process('CallProcess', 'process_import_bids', [], [], ...
    'bidsdir',          {BidsDir, 'BIDS'}, ...
    'selectsubj',       SubjectName, ...
    'downsamplemethod', 'icosphere', ...
    'icolevel',         IcoLevel, ...
    'nvertices',        15000, ...        % ignored on the icosphere path
    'channelalign',     0);

% process_import_bids returns [] -- or a PARTIAL list (rest only: 20 OMEGA controls, array
% 63280594, lost their task-noise run this way) -- even when every recording WAS imported (raw
% CTF links). So the handle list is ALWAYS rebuilt from the DB, scanning ALL of the subject's
% studies (NOT just @intra_subject, which never holds raw recordings).
nReturned = numel(sFilesRaw);
[sSubject, iSubject] = bst_get('Subject', SubjectName); %#ok<ASGLU>
if ~isempty(iSubject)
    sFilesDb = struct('FileName', {});
    sStudies = bst_get('StudyWithSubject', sSubject.FileName);   % all studies
    for iSt = 1:numel(sStudies)
        for iData = 1:numel(sStudies(iSt).Data)
            sFilesDb(end+1).FileName = sStudies(iSt).Data(iData).FileName; %#ok<AGROW>
        end
    end
    fprintf('Recordings: %d returned by process_import_bids, %d in the protocol DB.\n', nReturned, numel(sFilesDb));
    if numel(sFilesDb) >= nReturned
        sFilesRaw = sFilesDb;
    end
end
if isempty(sFilesRaw)
    error('BIDS import produced no data files for %s', SubjectName);
end

% Verify cortex exists (icosphere import must have read FreeSurfer derivatives)
sSubjects = bst_get('ProtocolSubjects');
for iSub = 1:numel(sSubjects.Subject)
    nm = sSubjects.Subject(iSub).Name; iC = sSubjects.Subject(iSub).iCortex;
    if (isempty(iC) || iC == 0) && ~strcmpi(nm,'@default_subject') && ~startsWith(nm,'sub-emptyroom')
        warning('No cortex for %s — check derivatives/freesurfer. Source stage will fail.', nm);
    end
end

% Head points cleanup + continuous CTF (resample happens in preprocess so the
% pre-resample raw handle survives for the 3-way folder cleanup there).
sFilesRaw = bst_process('CallProcess', 'process_headpoints_remove', sFilesRaw, [], 'zlimit', 0);
sFilesRaw = bst_process('CallProcess', 'process_headpoints_refine', sFilesRaw, []);
sFilesRaw = bst_process('CallProcess', 'process_ctf_convert',       sFilesRaw, [], 'rectype', 2);
end


%% ########################################################################
%% MODULE: PREPROCESS  (notch harmonics, hp 0.3, SSP, bad segments)
%% ########################################################################
function [sFilesRest, sFilesBand] = pa_preprocess(sFilesRaw)
% Resample 2400->1200 FIRST (lossless for <=300 Hz; halves samples/file size),
% then notch 60/120/180/240/300 + high-pass 0.3 Hz.
sFilesResample = bst_process('CallProcess', 'process_resample', sFilesRaw, [], 'freq', 1200);
sFilesNotch = bst_process('CallProcess', 'process_notch', sFilesResample, [], ...
    'freqlist', [60 120 180 240 300], 'sensortypes', 'MEG, EEG', 'read_all', 1);
sFilesBand  = bst_process('CallProcess', 'process_bandpass', sFilesNotch, [], ...
    'sensortypes', 'MEG, EEG', 'highpass', 0.3, 'lowpass', 0, ...
    'attenuation', 'strict', 'mirror', 0, 'useold', 0, 'read_all', 1);

% PSD (Welch) over the FULL recording, FULL spectrum (timewindow []=whole file,
% Freqs []=all bins). Sensor-level QC of the cleaned band-pass data.
bst_process('CallProcess', 'process_psd', sFilesBand, [], ...
    'timewindow', [], 'win_length', 4, 'win_overlap', 50, 'sensortypes', 'MEG, EEG', ...
    'edit', struct('Comment','Power', 'TimeBands',[], 'Freqs',[], ...
        'ClusterFuncTime','none', 'Measure','power', 'Output','all', 'SaveKernel',0));

% Drop intermediates (raw + resample + notch); keep band-pass + its PSD
bst_process('CallProcess', 'process_delete', [sFilesRaw, sFilesResample, sFilesNotch], [], 'target', 2);

% --- artifact cleaning on task-rest ---
sFilesRest = bst_process('CallProcess', 'process_select_tag', sFilesBand, [], ...
    'tag', 'task-rest', 'search', 1, 'select', 1);
bst_process('CallProcess', 'process_evt_detect_ecg', sFilesRest, [], ...
    'channelname', 'ECG',  'eventname', 'cardiac');
bst_process('CallProcess', 'process_evt_detect_eog', sFilesRest, [], ...
    'channelname', 'VEOG', 'eventname', 'blink');
bst_process('CallProcess', 'process_evt_detect_eog', sFilesRest, [], ...
    'channelname', 'HEOG', 'eventname', 'saccade');
bst_process('CallProcess', 'process_evt_remove_simult', sFilesRest, [], ...
    'remove', 'cardiac', 'target', 'blink', 'dt', 0.25, 'rename', 0);
bst_process('CallProcess', 'process_ssp_ecg', sFilesRest, [], ...
    'eventname', 'cardiac', 'sensortypes', 'MEG', 'usessp', 1, 'select', 1);
bst_process('CallProcess', 'process_ssp_eog', sFilesRest, [], ...
    'eventname', 'blink',   'sensortypes', 'MEG', 'usessp', 1, 'select', 1);
bst_process('CallProcess', 'process_ssp_eog', sFilesRest, [], ...
    'eventname', 'saccade', 'sensortypes', 'MEG', 'usessp', 1, 'select', 1);
% (No automatic bad-segment detection — bad segments are marked manually during
%  QC/validation in the GUI before any real analysis.)
end


%% ########################################################################
%% MODULE: SOURCE  (noise cov, OS head model, nxr manifold, dSPM + Dirac dSPM)
%% ########################################################################
function sSrc = pa_source(sFilesRest, sFilesBand)
% ⚠ bst_process NEVER throws: a failing process writes an 'error' line into the
% Brainstorm report and returns. Every call below goes through pa_call, which turns
% such an entry into a MATLAB error, and every product (noise cov, head model, dSPM
% kernel, Dirac kernel) is asserted to exist for every rest run. A subject without
% a head model or a kernel is a FAILED task, never a completed one. (OMEGA bst-meg
% 2026-10-06: 104 of 524 archives were promoted with no head model / no kernel.)
if isempty(sFilesRest)
    error('pa_source: no task-rest recording to source-localise.');
end
sFilesNoise = bst_process('CallProcess', 'process_select_tag', sFilesBand, [], ...
    'tag', 'task-noise', 'search', 1, 'select', 1);
if isempty(sFilesNoise)
    error('No task-noise recording found for noise covariance.');
end

% --- Noise covariance: one per noise recording, then assigned explicitly ---
% ⚠ process_noisecov 'copymatch' copies a noise cov only to folders where that noise
% recording is STRICTLY closer in DateOfStudy than every other noise recording. Two
% noise runs on the same day tie, so NOTHING is copied and the inverse has no noise
% cov (OMEGA PD/MNI subjects with task-noise run-01 + run-02). So: compute without
% copying, then give each rest run the noise of its own session (same ses- token),
% else the closest date, ties to the first.
for iN = 1:numel(sFilesNoise)
    pa_call('process_noisecov', sFilesNoise(iN), ...
        'target', 1, 'dcoffset', 1, 'identity', 0, ...
        'copycond', 0, 'copysubj', 0, 'copymatch', 0, 'replacefile', 1);
end
for iR = 1:numel(sFilesRest)
    iN = pa_match_noise(sFilesRest(iR).FileName, {sFilesNoise.FileName});
    [~, iSrcCh] = bst_get('ChannelForStudy', pa_istudy(sFilesNoise(iN).FileName));
    [~, iDstCh] = bst_get('ChannelForStudy', pa_istudy(sFilesRest(iR).FileName));
    if ~isequal(iSrcCh, iDstCh)
        db_set_noisecov(iSrcCh, iDstCh, 0, 1);
    end
    fprintf('Noise cov: %s <- %s\n', sFilesRest(iR).FileName, sFilesNoise(iN).FileName);
    pa_sync_badchannels(sFilesRest(iR).FileName);
end
pa_assert_study(sFilesRest, 'NoiseCov', 'noise covariance');

% ONE overlapping-spheres head model feeds BOTH inverses. MEG ONLY: the EEG default
% of process_headmodel is OpenMEEG BEM (eeg=3), which needs an inner skull. ⚠ OMEGA
% recordings that carry EEG-typed channels therefore failed the WHOLE head model
% with "OpenMEEG: Inner skull surface not available" — silently.
% Registration fallback: when the head-point refinement (ICP, pa_import) pulls the helmet into
% the head ("sensors ... inside the brain volume"; OMEGA sub-0153, sub-CONP0078 = PREVENT-AD
% sub-MTL0160, sub-MNI0079), remove that refinement -- back to the head-coil/fiducial alignment
% -- and retry ONCE. A second failure is a data problem (MRI or digitisation) and fails the task.
try
    pa_call('process_headmodel', sFilesRest, ...
        'sourcespace', 1, 'meg', 3, 'eeg', 1, 'ecog', 1, 'seeg', 1);
catch err
    if isempty(strfind(err.message, 'inside the brain')), rethrow(err); end
    fprintf('REGISTRATION FALLBACK: %s\n  -> removing the head-point refinement, retrying the head model (fiducials only)\n', err.message);
    pa_call('process_adjust_coordinates', sFilesRest, ...
        'reset', 0, 'head', 0, 'points', 1, 'remove', 1, 'display', 0);
    pa_call('process_headmodel', sFilesRest, ...
        'sourcespace', 1, 'meg', 3, 'eeg', 1, 'ecog', 1, 'seeg', 1);
    fprintf('REGISTRATION FALLBACK: head model OK without the head-point refinement\n');
end
pa_assert_study(sFilesRest, 'HeadModel', 'head model');

% nxr manifold backbone (find-or-create) — needs nxr-compute plugin
pa_call('process_tess_manifold', sFilesRest, ...
    'gauge', 'trivial', 'forcerecompute', 0);

% Standard dSPM (unconstrained, kernel only)
sSrc = pa_call('process_inverse_2018', sFilesRest, ...
    'output', 2, ...
    'inverse', struct('Comment','dSPM: MEG', 'InverseMethod','minnorm', ...
        'InverseMeasure','dspm2018', 'SourceOrient',{{'free'}}, 'Loose',0.2, ...
        'UseDepth',1, 'WeightExp',0.5, 'WeightLimit',10, 'NoiseMethod','reg', ...
        'NoiseReg',0.1, 'SnrMethod','fixed', 'SnrRms',1e-06, 'SnrFixed',3, ...
        'ComputeKernel',1, 'DataTypes',{{'MEG'}}));

% Custom Dirac dSPM (Dirac eigenbasis; same OS head model) — needs nxr-compute
pa_call('process_inverse_dirac', sFilesRest, ...
    'measure','dspm2018', 'snr',3, 'nmodes',400, 'tau',0.5, ...
    'noisereg',0.1, 'sensortypes','MEG');

pa_assert_kernels(sFilesRest);
end


%% ===== SOURCE helper: bad channels of the noise run are bad in the rest run =====
function pa_sync_badchannels(restFile)
% A channel that was bad in the noise recording has an all-zero row/column in the noise
% covariance; if it is good in the rest run, process_inverse_2018 refuses with "Bad channels in
% noise covariance are different from bad channels in recordings" (OMEGA sub-PD0457). Mark such
% MEG channels bad in the rest run too: they cannot be whitened, so they are left out.
[~, iStudy] = bst_get('AnyFile', restFile);
[sCh, iChStudy] = bst_get('ChannelForStudy', iStudy);
sStCh = bst_get('Study', iChStudy);
if isempty(sStCh.NoiseCov), return; end
N = load(file_fullpath(sStCh.NoiseCov(1).FileName), 'NoiseCov');
if isempty(N.NoiseCov), return; end
ChannelMat = in_bst_channel(sCh.FileName);
D = in_bst_data(restFile, 'ChannelFlag');
iGood = good_channel(ChannelMat.Channel, D.ChannelFlag, 'MEG');
iBad = intersect(find(~any(N.NoiseCov, 1) & ~any(N.NoiseCov, 2)'), iGood);
if isempty(iBad), return; end
names = strjoin({ChannelMat.Channel(iBad).Name}, ', ');
fprintf('Bad in the noise run, now bad in %s: %s\n', restFile, names);
pa_call('process_channel_setbad', restFile, 'sensortypes', names);
end


%% ===== SOURCE helpers: make bst_process fail loudly =====
function sOut = pa_call(procName, sInputs, varargin)
% Run one process; a new 'error' entry in the Brainstorm report becomes a MATLAB error.
global GlobalData
n0 = 0;
if ~isempty(GlobalData) && iscell(GlobalData.ProcessReports.Reports)
    n0 = size(GlobalData.ProcessReports.Reports, 1);
end
sOut = bst_process('CallProcess', procName, sInputs, [], varargin{:});
R = GlobalData.ProcessReports.Reports;
if iscell(R) && size(R,1) > n0
    isErr = strcmpi(R(n0+1:end,1), 'error');
    if any(isErr)
        msgs = R(n0 + find(isErr), 4);
        msgs = cellfun(@(m) strtok(char(m), char(10)), msgs, 'UniformOutput', false);
        error('%s failed: %s', procName, strjoin(unique(msgs, 'stable'), ' | '));
    end
end
end

function iStudy = pa_istudy(FileName)
% Current study index of a file (indices in sFiles go stale after process_delete).
[~, iStudy] = bst_get('AnyFile', FileName);
if isempty(iStudy); error('File not in the database: %s', FileName); end
end

function iN = pa_match_noise(restFile, noiseFiles)
% Noise recording for one rest run: same ses- token first, then closest DateOfStudy.
ses = regexp(restFile, 'ses-[A-Za-z0-9]+', 'match', 'once');
cand = 1:numel(noiseFiles);
if ~isempty(ses)
    same = find(contains(noiseFiles, ses));
    if ~isempty(same); cand = same; end
end
iN = cand(1);
if numel(cand) > 1
    dR = pa_date(restFile);
    if ~isempty(dR)
        d = arrayfun(@(k) pa_date(noiseFiles{k}), cand, 'UniformOutput', false);
        ok = ~cellfun(@isempty, d);
        if any(ok)
            dd = inf(size(cand)); dd(ok) = abs([d{ok}] - dR);
            [~, k] = min(dd);    % ties -> first (lowest run)
            iN = cand(k);
        end
    end
end
end

function d = pa_date(FileName)
sStudy = bst_get('AnyFile', FileName);
d = [];
if ~isempty(sStudy) && ~isempty(sStudy.DateOfStudy)
    try, d = datenum(sStudy.DateOfStudy); catch, d = []; end %#ok<DATNM>
end
end

function pa_assert_study(sFiles, field, what)
% Every input's channel study must hold a non-empty <field> (NoiseCov | HeadModel).
for i = 1:numel(sFiles)
    iSt = pa_istudy(sFiles(i).FileName);
    [~, iChSt] = bst_get('ChannelForStudy', iSt);   % noise cov + head model live here
    sStudy = bst_get('Study', iChSt);
    v = sStudy.(field);
    if isempty(v) || isempty(v(1).FileName)
        error('No %s for %s — source step failed.', what, sFiles(i).FileName);
    end
end
end

function pa_assert_kernels(sFiles)
% Every rest run must have a standard dSPM kernel AND a Dirac kernel.
for i = 1:numel(sFiles)
    iSt = pa_istudy(sFiles(i).FileName);
    [~, iChSt] = bst_get('ChannelForStudy', iSt);   % shared kernels live with the channel file
    fn = {};
    for jSt = unique([iSt, iChSt])
        sStudy = bst_get('Study', jSt);
        if ~isempty(sStudy.Result); fn = [fn, {sStudy.Result.FileName}]; end %#ok<AGROW>
    end
    isK  = contains(fn, 'KERNEL');
    isDi = contains(fn, 'Dirac', 'IgnoreCase', true);
    if ~any(isK & ~isDi)
        error('No dSPM kernel for %s — source step failed.', sFiles(i).FileName);
    end
    if ~any(isDi)
        error('No Dirac dSPM result for %s — source step failed.', sFiles(i).FileName);
    end
    fprintf('Kernels OK: %s (%d results)\n', sFiles(i).FileName, numel(fn));
end
end


%% ########################################################################
%% MODULE: TIMEFREQ  (band-power maps on dSPM sources, projected to template)
%% ########################################################################
function pa_power(sFilesRest)
% Resolve the standard-dSPM source links for these rest recordings.
srcLinks = {};
for iR = 1:numel(sFilesRest)
    sStudyR = bst_get('AnyFile', sFilesRest(iR).FileName);
    if isempty(sStudyR) || isempty(sStudyR.Result); continue; end
    for iRes = 1:numel(sStudyR.Result)
        fn = sStudyR.Result(iRes).FileName;
        if ~isempty(strfind(fn,'link|')) && isempty(strfind(fn,'DiracEig'))
            srcLinks{end+1} = fn; %#ok<AGROW>
        end
    end
end
if isempty(srcLinks)
    error('No dSPM source links found for the rest recordings; cannot compute power maps.');
end
sP = bst_process('CallProcess', 'process_psd', srcLinks, [], ...
    'timewindow', [0 100], 'win_length', 4, 'win_overlap', 50, 'scoutfunc', 1, ...
    'edit', struct('Comment','Power,FreqBands', 'TimeBands',[], ...
        'Freqs',{{'delta','2, 4','mean'; 'theta','5, 7','mean'; 'alpha','8, 12','mean'; ...
                  'beta','15, 29','mean'; 'gamma1','30, 59','mean'; 'gamma2','60, 90','mean'}}, ...
        'ClusterFuncTime','none', 'Measure','power', 'Output','all', 'SaveKernel',0));
sP = bst_process('CallProcess', 'process_tf_norm', sP, [], 'normalize','relative', 'overwrite',0);
sP = bst_process('CallProcess', 'process_project_sources', sP, [], 'headmodeltype','surface');
sP = bst_process('CallProcess', 'process_ssmooth_surfstat', sP, [], 'fwhm',3, 'overwrite',1);
bst_process('CallProcess', 'process_average', sP, [], ...
    'avgtype',1, 'avg_func',1, 'weighted',0, 'matchrows',0, 'iszerobad',0);
end


%% ===== helper: expand leading ~ in a path (Java IO treats ~ as literal) =====
function pth = expand_tilde(pth)
if ~isempty(pth) && pth(1) == '~'
    home = char(java.lang.System.getProperty('user.home'));
    if numel(pth) == 1
        pth = home;
    elseif pth(2) == filesep
        pth = fullfile(home, pth(3:end));
    end
end
end
