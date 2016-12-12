function tutorial_visual_single(tutorial_dir, reports_dir, iSubjStart)
% TUTORIAL_VISUAL_SINGLE: Runs the Brainstorm/SPM group analysis pipeline (single subject).
%
% ONLINE TUTORIALS: http://neuroimage.usc.edu/brainstorm/Tutorials/VisualSingle
%
% INPUTS:
%    - tutorial_dir: Directory containing the folder sample_group.
%       |- sample_group
%           |- anatomy/subXXX/             : Segmentation folders generated with FreeSurfer (downloaded from Brainstorm website: sample_group_freesurfer.zip)
%           |- ds117/subXXX/MEG/*_sss.fif  : MEG+EEG recordings downloaded from https://openfmri.org/dataset/ds000117/
%           |- emptyroom/090707_raw_st.fif : Empty room measurements (http://openfmri.s3.amazonaws.com/tarballs/ds117_metadata.tgz)
%    - reports_dir: If defined, exports all the reports as HTML to this folder
%    - iSubjStart : Index of the first subject to process
%           => If the script crashes, you can re-run it starting from the last subject

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Author: Francois Tadel, Elizabeth Bock, 2016


%% ===== SCRIPT VARIABLES =====
% Full list of subjects to process
SubjectNames = {'sub001', 'sub002', 'sub003', 'sub004', 'sub005', 'sub006', 'sub007', 'sub008', 'sub009', 'sub010', ...
                'sub011', 'sub012', 'sub013', 'sub014', 'sub015', 'sub016', 'sub017', 'sub018', 'sub019'};
SubjectNoise = 'emptyroom';
% Bad channels {iSubj} = {Run01, Run02, Run03, Run04, Run05, Run06}
BadChannels{1}   = {[], {'EEG009','EEG039'}, [], [], [], []};
BadChannels{2}   = {'EEG016', 'EEG070', 'EEG050',{'EEG008','EEG050'}, [], []};
BadChannels{3}   = {{'EEG027', 'EEG030', 'EEG038'}, 'EEG010', 'EEG010', 'EEG010', 'EEG010', 'EEG010'};
BadChannels{4}   = {{'EEG008','EEG017'}, {'EEG008','EEG017'}, {'EEG008','EEG017'}, {'EEG008','EEG017'}, {'EEG008','EEG017','EEG001'}, {'EEG008','EEG017','EEG020'}};
BadChannels{5}   = {[], [], [], [], [], []};
BadChannels{6}   = {[], 'EEG004', 'EEG004', [], 'EEG004', 'EEG004'};
BadChannels{7}   = {{'EEG038'}, {'EEG038','EEG001','EEG016'}, {'EEG038','EEG001','EEG016'}, {'EEG038','EEG001'}, {'EEG038','EEG001','EEG016'}, {'EEG038','EEG001','EEG016'}};
BadChannels{8}   = {'EEG001', 'EEG001', [], [], [], []};
BadChannels{9}   = {'EEG068', [], 'EEG004', [], [], []};
BadChannels{10}  = {[], [], {'EEG004','EEG008'}, {'EEG004','EEG008'},{'EEG004','EEG008','EEG043','EEG045','EEG047'}, {'EEG004','EEG008'}};
BadChannels{11}  = {[], [], [], [], [], []};
BadChannels{12}  = {[], [], [], [], [], []};
BadChannels{13}  = {{'EEG010','EEG050'}, 'EEG050', 'EEG050', 'EEG050', 'EEG050', 'EEG050'};
BadChannels{14}  = {{'EEG024','EEG057'}, {'EEG024','EEG057'}, {'EEG024','EEG057'}, {'EEG024','EEG057','EEG070'}, {'EEG024','EEG057'}, {'EEG024','EEG057','EEG070'}};
BadChannels{15}  = {'EEG009', 'EEG009', {'EEG009','EEG057','EEG69'}, 'EEG009', {'EEG009','EEG044'}, {'EEG009','EEG044'}};
BadChannels{16}  = {[], [], [], [], [], []};
BadChannels{17}  = {'EEG029', 'EEG029', 'EEG029', {'EEG004','EEG008','EEG016','EEG029'}, {'EEG004','EEG008','EEG016','EEG029'}, {'EEG004','EEG008','EEG016','EEG029'}};
BadChannels{18}  = {'EEG038', 'EEG038', 'EEG038', 'EEG038', {'EEG054','EEG038'}, 'EEG038'};
BadChannels{19}  = {'EEG008', 'EEG008', 'EEG008', 'EEG008', 'EEG008', 'EEG008'};
% SSP components to remove {iSubj} = {sRun01, sRun02, sRun03, sRun03, sRun04, sRun05, sRun06},   sRun0X={ECG_GRAD,ECG_MAG}
SspSelect{1}  = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{2}  = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}}; 
SspSelect{3}  = {{1,1},         {1,1},   {1,1},   {1,1},   {3,1},   {1,1}};
SspSelect{4}  = {{[],1},        {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{5}  = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{6}  = {{1,1},         {1,1},   {[],1},  {[],1},  {1,1},   {1,1}};
SspSelect{7}  = {{[],1},        {[],1},  {[],1},  {[],1},  {[],1},  {[],1}};
SspSelect{8}  = {{2,1},         {1,1},   {1,1},   {[],1},  {1,1},   {1,1}};
SspSelect{9}  = {{2,1},         {2,1},   {1,1},   {2,1},   {1,1},   {2,1}};
SspSelect{10} = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{11} = {{1,1},         {1,1},   {[],1},  {2,1},   {1,1},   {2,1}};
SspSelect{12} = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{13} = {{[],1},        {[],1},  {[],1},  {[],1},  {[],[]}, {[],[]}};
SspSelect{14} = {{[1,2],[1,2]}, {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{15} = {{[],[]},       {[],[]}, {[],[]}, {[],[]}, {[],[]}, {[],[]}};
SspSelect{16} = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{17} = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{18} = {{1,1},         {1,1},   {1,1},   {1,1},   {1,1},   {1,1}};
SspSelect{19} = {{1,1},         {1,1},   {1,1},   {2,1},   {1,1},   {1,1}};


%% ===== CREATE PROTOCOL =====
% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end
% First subject to process
if (nargin < 3) || isempty(iSubjStart)
    iSubjStart = 1;
end
% Output folder for reports
if (nargin < 2) || isempty(reports_dir) || ~isdir(reports_dir)
    reports_dir = [];
end
% You have to specify the folder in which the tutorial dataset is unzipped
if (nargin < 1) || isempty(tutorial_dir) || ~file_exist(tutorial_dir)
    error('The first argument must be the full path to the tutorial folder.');
end
% The protocol name has to be a valid folder name (no spaces, no weird characters...)
ProtocolName = 'TutorialVisual';
% If starting from the first subject: delete the protocol
if (iSubjStart == 1)
    % Delete existing protocol
    gui_brainstorm('DeleteProtocol', ProtocolName);
    % Create new protocol
    gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);
end
% Set visualization filters: 40Hz low-pass, no high-pass
panel_filter('SetFilters', 1, 40, 0, [], 0, [], 0, 0);


%% ===== EMPTY ROOM RECORDINGS =====
% If starting from the first subject: import the noise covariance
if (iSubjStart == 1)
    % Use the same noise recordings for all the files
    NoiseFile  = fullfile(tutorial_dir, 'sample_group', 'emptyroom', '090707_raw_st.fif');
    % Process: Create link to raw file
    sFileNoise = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
        'subjectname',    SubjectNoise, ...
        'datafile',       {NoiseFile, 'FIF'}, ...
        'channelreplace', 1, ...
        'channelalign',   0);
    % Process: Notch filter: 50Hz 100Hz 150Hz 200Hz
    sFileNoiseClean = bst_process('CallProcess', 'process_notch', sFileNoise, [], ...
            'freqlist',    [50, 100, 150, 200], ...
            'sensortypes', 'MEG, EEG', ...
            'read_all',    0);
    % Process: Compute noise covariance
    bst_process('CallProcess', 'process_noisecov', sFileNoiseClean, [], ...
        'baseline',    [], ...
        'sensortypes', 'MEG, EEG', ...
        'target',      1, ...  % Noise covariance     (covariance over baseline time window)
        'dcoffset',    1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',    0, ...
        'copycond',    0, ...
        'copysubj',    0, ...
        'replacefile', 1);  % Replace
% If the noise recordings are already imported: just select them from the database
else
    % Process: Select data files in: emptyroom/*
    sFileNoiseClean = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   SubjectNoise);
    % Process: Select file comments with tag: Avg
    sFileNoiseClean = bst_process('CallProcess', 'process_select_tag', sFileNoiseClean, [], ...
        'tag',    'notch');  % Select only the files with the tag
    % Check if it was already processed
    if isempty(sFileNoiseClean)
        error('Noise recordings have not been imported yet. Re-run the script from subject #1.');
    end
    % Select only the first file if there are multiple files
    sFileNoiseClean = sFileNoiseClean(1);
end


%% ===== PRE-PROCESS AND IMPORT =====
for iSubj = iSubjStart:length(SubjectNames)
    % Start a new report (one report per subject)
    bst_report('Start');
    
    % If subject already exists: delete it
    [sSubject, iSubject] = bst_get('Subject', SubjectNames{iSubj});
    if ~isempty(sSubject)
        db_delete_subjects(iSubject);
    end
    
    % ===== FILES TO IMPORT =====
    % Build the path of the files to import
    AnatDir    = fullfile(tutorial_dir, 'sample_group', 'freesurfer', SubjectNames{iSubj});
    DataDir    = fullfile(tutorial_dir, 'sample_group', 'ds117', SubjectNames{iSubj}, 'MEG');
    % Check if the folder contains the required files
    if ~file_exist(AnatDir)
        error(['The folder "' AnatDir '" does not exist.']);
    end
    if ~file_exist(DataDir)
        error(['The folder "' DataDir '" does not exist.']);
    end

    % ===== ANATOMY =====
    % Process: Import anatomy folder
    bst_process('CallProcess', 'process_import_anatomy', [], [], ...
        'subjectname', SubjectNames{iSubj}, ...
        'mrifile',     {AnatDir, 'FreeSurfer'}, ...
        'nvertices',   15000);
    
    % ===== PROCESS EACH RUN =====
    for iRun = 1:6
        % Files to import
        FifFile = bst_fullfile(DataDir, sprintf('run_%02d_sss.fif', iRun));
        BadFile = fullfile(tutorial_dir, 'sample_group', 'brainstorm', 'bad_segments', sprintf('sub%03d_run_%02d_sss_notch_events.mat',iSubj,iRun));

        % ===== LINK CONTINUOUS FILE =====
        % Process: Create link to raw file
        sFileRaw = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
            'subjectname',    SubjectNames{iSubj}, ...
            'datafile',       {FifFile, 'FIF'}, ...
            'channelreplace', 1, ...
            'channelalign',   0);

        % ===== PREPARE CHANNEL FILE =====
        % Process: Set channels type
        bst_process('CallProcess', 'process_channel_settype', sFileRaw, [], ...
            'sensortypes', 'EEG061, EEG064', ...
            'newtype',     'NOSIG');
        bst_process('CallProcess', 'process_channel_settype', sFileRaw, [], ...
            'sensortypes', 'EEG062', ...
            'newtype',     'EOG');
        bst_process('CallProcess', 'process_channel_settype', sFileRaw, [], ...
            'sensortypes', 'EEG063', ...
            'newtype',     'ECG');

        % Process: Remove head points
        sFileRaw = bst_process('CallProcess', 'process_headpoints_remove', sFileRaw, [], ...
            'zlimit', 0);
        % Process: Refine registration
        sFileRaw = bst_process('CallProcess', 'process_headpoints_refine', sFileRaw, []);
        % Process: Project electrodes on scalp
        sFileRaw = bst_process('CallProcess', 'process_channel_project', sFileRaw, []);

        % Process: Snapshot: Sensors/MRI registration
        bst_process('CallProcess', 'process_snapshot', sFileRaw, [], ...
            'target',   1, ...  % Sensors/MRI registration
            'modality', 1, ...  % MEG (All)
            'orient',   1, ...  % left
            'comment',  sprintf('MEG/MRI Registration: Subject #%d, Run #%d', iSubj, iRun));
        bst_process('CallProcess', 'process_snapshot', sFileRaw, [], ...
            'target',   1, ...  % Sensors/MRI registration
            'modality', 4, ...  % EEG
            'orient',   1, ...  % left
            'comment',  sprintf('EEG/MRI Registration: Subject #%d, Run #%d', iSubj, iRun));

        % ===== IMPORT TRIGGERS =====
        % Process: Read from channel
        bst_process('CallProcess', 'process_evt_read', sFileRaw, [], ...
            'stimchan',  'STI101', ...
            'trackmode', 2, ...  % Bit: detect the changes for each bit independently
            'zero',      0);
        % Process: Group by name
        bst_process('CallProcess', 'process_evt_groupname', sFileRaw, [], ...
            'combine', 'Unfamiliar=3,4', ...
            'dt',      0, ...
            'delete',  1);
        % Process: Rename event
        bst_process('CallProcess', 'process_evt_rename', sFileRaw, [], ...
            'src',  '3', ...
            'dest', 'Famous');
        % Process: Rename event
        bst_process('CallProcess', 'process_evt_rename', sFileRaw, [], ...
            'src',  '5', ...
            'dest', 'Scrambled');
        % Process: Add time offset
        bst_process('CallProcess', 'process_evt_timeoffset', sFileRaw, [], ...
            'info',      [], ...
            'eventname', 'Famous, Unfamiliar, Scrambled', ...
            'offset',    0.0345);
        % Process: Delete events
        bst_process('CallProcess', 'process_evt_delete', sFileRaw, [], ...
            'eventname', '1,2,6,7,8,9,10,11,12,13,14,15,16');
        % Process: Detect cHPI activity (Elekta):STI201
        bst_process('CallProcess', 'process_evt_detect_chpi', sFileRaw, [], ...
            'eventname',   'chpi_bad', ...
            'channelname', 'STI201', ...
            'method',      'off');  % Mark as bad when the HPI coils are OFF
        
        % ===== FREQUENCY FILTERS =====
        % Process: Notch filter: 50Hz 100Hz 150Hz 200Hz
        sFileClean = bst_process('CallProcess', 'process_notch', sFileRaw, [], ...
            'freqlist',    [50, 100, 150, 200], ...
            'sensortypes', 'MEG, EEG', ...
            'read_all',    0);
        % Process: Power spectrum density (Welch)
        sFilesPsd = bst_process('CallProcess', 'process_psd', [sFileRaw, sFileClean], [], ...
            'timewindow',  [], ...
            'win_length',  4, ...
            'win_overlap', 50, ...
            'sensortypes', 'MEG, EEG', ...
            'edit', struct(...
                 'Comment',         'Power', ...
                 'TimeBands',       [], ...
                 'Freqs',           [], ...
                 'ClusterFuncTime', 'none', ...
                 'Measure',         'power', ...
                 'Output',          'all', ...
                 'SaveKernel',      0));
        % Process: Snapshot: Frequency spectrum
        bst_process('CallProcess', 'process_snapshot', sFilesPsd, [], ...
            'target',   10, ...  % Frequency spectrum
            'comment',  sprintf('Power spctrum: Subject #%d, Run #%d', iSubj, iRun));

        % ===== BAD CHANNELS =====
        if ~isempty(BadChannels{iSubj}{iRun})
            % Process: Set bad channels
            bst_process('CallProcess', 'process_channel_setbad', sFileClean, [], ...
                'sensortypes', BadChannels{iSubj}{iRun});
        end

        % ===== EEG REFERENCE =====
        % Process: Re-reference EEG
        bst_process('CallProcess', 'process_eegref', sFileClean, [], ...
            'eegref',      'AVERAGE', ...
            'sensortypes', 'EEG');
    
        % ===== DETECT ARTIFACTS ======
        % Process: Detect heartbeats
        bst_process('CallProcess', 'process_evt_detect_ecg', sFileClean, [], ...
            'channelname', 'EEG063', ...
            'timewindow',  [], ...
            'eventname',   'cardiac');
        % Different amplitude thresholds for different subjects
        if strcmpi(SubjectNames{iSubj}, 'sub006')
            thresholdMAX = 50;
        else
            thresholdMAX = 100;
        end
        % Process: Detect: blink_BAD    -   Detects all events where the amplitude exceeds 100uV
        bst_process('CallProcess', 'process_evt_detect_threshold', sFileClean, [], ...
            'eventname',    'blink_BAD', ...
            'channelname',  'EEG062', ...
            'timewindow',   [], ...
            'thresholdMAX', thresholdMAX, ...
            'units',        3, ...  % uV (10^-6)
            'bandpass',     [0.3, 20], ...
            'isAbsolute',   1, ...
            'isDCremove',   0);

        % ===== SSP COMPUTATION =====
        % Process: SSP ECG: cardiac
        bst_process('CallProcess', 'process_ssp_ecg', sFileClean, [], ...
            'eventname',   'cardiac', ...
            'sensortypes', 'MEG GRAD', ...
            'usessp',      1, ...
            'select',      SspSelect{iSubj}{iRun}{1}); 
        bst_process('CallProcess', 'process_ssp_ecg', sFileClean, [], ...
            'eventname',   'cardiac', ...
            'sensortypes', 'MEG MAG', ...
            'usessp',      1, ...
            'select',      SspSelect{iSubj}{iRun}{2});
        % Process: Snapshot: SSP projectors
        bst_process('CallProcess', 'process_snapshot', sFileClean, [], ...
            'target',   2, ...
            'comment',  sprintf('Subject #%d, Run #%d', iSubj, iRun));   % SSP projectors
        
        % ===== IMPORT BAD EVENTS =====
        % Process: Import from file
        bst_process('CallProcess', 'process_evt_import', sFileClean, [], ...
            'evtfile', {BadFile, 'BST'});
        
        % ===== IMPORT TRIALS =====
        % Process: Import MEG/EEG: Events
        sFilesEpochs = bst_process('CallProcess', 'process_import_data_event', sFileClean, [], ...
            'subjectname', SubjectNames{iSubj}, ...
            'condition',   '', ...
            'eventname',   'Famous, Scrambled, Unfamiliar', ...
            'timewindow',  [], ...
            'epochtime',   [-0.5, 1.2], ...
            'createcond',  0, ...
            'ignoreshort', 1, ...
            'usectfcomp',  1, ...
            'usessp',      1, ...
            'freq',        [], ...
            'baseline',    [-0.5, -0.0009]);
        
        % ===== AVERAGE: RUN =====
        % Process: Average: By trial group (folder average)
        sFilesAvg = bst_process('CallProcess', 'process_average', sFilesEpochs, [], ...
            'avgtype',    5, ...  % By trial group (folder average)
            'avg_func',   1, ...  % Arithmetic average:  mean(x)
            'weighted',   0, ...
            'keepevents', 0);
        % Process: Snapshot: Recordings time series
        bst_process('CallProcess', 'process_snapshot', sFilesAvg, [], ...
            'target',   5, ...  % Recordings time series
            'modality', 4, ...  % EEG
            'time',     0.11, ...
            'Comment',  sprintf('Subject #%d, Run #%d', iSubj, iRun));
        % Process: Snapshot: Recordings time series
        bst_process('CallProcess', 'process_snapshot', sFilesAvg, [], ...
            'target',   6, ...  % Recordings topography (one time)
            'modality', 4, ...  % EEG
            'time',     0.11, ...
            'Comment',  sprintf('Subject #%d, Run #%d', iSubj, iRun));

        % ===== COPY NOISECOV: MEG =====
        % Copy noise covariance to current study
        isDataCov = 0;
        AutoReplace = 1;
        db_set_noisecov(sFileNoiseClean.iStudy, sFilesAvg(1).iStudy, isDataCov, AutoReplace);
        
        % ===== COMPUTE NOISECOV: EEG =====
        % Process: Compute covariance (noise or data)
        bst_process('CallProcess', 'process_noisecov', sFilesEpochs, [], ...
            'baseline',       [-0.5, -0.0009], ...
            'sensortypes',    'EEG', ...
            'target',         1, ...  % Noise covariance     (covariance over baseline time window)
            'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
            'identity',       0, ...
            'copycond',       0, ...
            'copysubj',       0, ...
            'replacefile',    2);  % Merge
    end
    
    % Save report
    ReportFile = bst_report('Save', []);
    if ~isempty(reports_dir) && ~isempty(ReportFile)
        bst_report('Export', ReportFile, bst_fullfile(reports_dir, ['report_' ProtocolName '_' SubjectNames{iSubj} '.html']));
    end
end



%% ===== SOURCE ESTIMATION =====
% Start a new report (one report for the source estimation of all the subjects)
bst_report('Start');
% Loop on the subjects: This loop is separated from the previous one, because we should 
% compute the BEM surfaces after importing all the runs, so that the registration is done 
% using the high resolution head surface, instead of the smooth scalp BEM layer.
for iSubj = 1:length(SubjectNames)
    % ===== BEM SURFACES =====
    % Process: Generate BEM surfaces
    bst_process('CallProcess', 'process_generate_bem', [], [], ...
        'subjectname', SubjectNames{iSubj}, ...
        'nscalp',      1082, ...
        'nouter',      642, ...
        'ninner',      642, ...
        'thickness',   4);

    % ===== SELECT ALL AVERAGES =====
    % Process: Select data files in: */*
    sFilesAvg = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   SubjectNames{iSubj});
    % Process: Select file comments with tag: Avg
    sFilesAvg = bst_process('CallProcess', 'process_select_tag', sFilesAvg, [], ...
        'tag',    'Avg');  % Select only the files with the tag

    % ===== COMPUTE HEAD MODELS =====
    % Process: Compute head model (only for the first run of the subject)
    bst_process('CallProcess', 'process_headmodel', sFilesAvg(1), [], ...
        'sourcespace', 1, ...  % Cortex surface
        'meg',         3, ...  % Overlapping spheres
        'eeg',         3, ...  % OpenMEEG BEM
        'ecog',        1, ...  % 
        'seeg',        1, ...  % 
        'openmeeg',    struct(...
             'BemSelect',    [1, 1, 1], ...
             'BemCond',      [1, 0.0125, 1], ...
             'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
             'BemFiles',     {{}}, ...
             'isAdjoint',    0, ...
             'isAdaptative', 1, ...
             'isSplit',      0, ...
             'SplitLength',  4000));
    % Get all the runs for this subject (ie the list of the study indices)
    iStudyOther = setdiff(unique([sFilesAvg.iStudy]), sFilesAvg(1).iStudy);
    % Copy the forward model file to the other runs
    sHeadmodel = bst_get('HeadModelForStudy', sFilesAvg(1).iStudy);
    for iStudy = iStudyOther
        db_add(iStudy, sHeadmodel.FileName);
    end
    
    % ===== COMPUTE SOURCES: MEG =====
    % Process: Compute sources [2016]
    sAvgSrcMeg = bst_process('CallProcess', 'process_inverse_2016', sFilesAvg, [], ...
        'output',  1, ...  % Kernel only: shared
        'inverse', struct(...
             'Comment',        'MN: MEG ALL', ...
             'InverseMethod',  'minnorm', ...
             'InverseMeasure', 'amplitude', ...
             'SourceOrient',   {{'fixed'}}, ...
             'Loose',          0.2, ...
             'UseDepth',       1, ...
             'WeightExp',      0.5, ...
             'WeightLimit',    10, ...
             'NoiseMethod',    'reg', ...
             'NoiseReg',       0.1, ...
             'SnrMethod',      'fixed', ...
             'SnrRms',         1e-06, ...
             'SnrFixed',       3, ...
             'ComputeKernel',  1, ...
             'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));
    % Process: Snapshot: Sources (one time)   -   Loop only to get a correct comment for the report
    for i = 1:length(sAvgSrcMeg)
        bst_process('CallProcess', 'process_snapshot', sAvgSrcMeg(i), [], ...
            'target',         8, ...  % Sources (one time)
            'orient',         4, ...  % bottom
            'time',           0.11, ...
            'threshold',      20, ...
            'Comment',        ['MEG sources: ' sFilesAvg(i).FileName]);
    end
    
    % ===== COMPUTE SOURCES: EEG =====
    % Process: Compute sources [2016]
    sAvgSrcEeg = bst_process('CallProcess', 'process_inverse_2016', sFilesAvg, [], ...
        'output',  1, ...  % Kernel only: shared
        'inverse', struct(...
             'Comment',        'MN: EEG', ...
             'InverseMethod',  'minnorm', ...
             'InverseMeasure', 'amplitude', ...
             'SourceOrient',   {{'fixed'}}, ...
             'Loose',          0.2, ...
             'UseDepth',       1, ...
             'WeightExp',      0.5, ...
             'WeightLimit',    10, ...
             'NoiseMethod',    'reg', ...
             'NoiseReg',       0.1, ...
             'SnrMethod',      'fixed', ...
             'SnrRms',         1e-06, ...
             'SnrFixed',       3, ...
             'ComputeKernel',  1, ...
             'DataTypes',      {{'EEG'}}));
    % Process: Snapshot: Sources (one time)   -   Loop only to get a correct comment for the report
    for i = 1:length(sAvgSrcEeg)
        bst_process('CallProcess', 'process_snapshot', sAvgSrcEeg(i), [], ...
            'target',         8, ...  % Sources (one time)
            'orient',         4, ...  % bottom
            'time',           0.11, ...
            'threshold',      10, ...
            'Comment',        ['EEG sources: ' sFilesAvg(i).FileName]);
    end
end
% Save report
ReportFile = bst_report('Save', []);
if ~isempty(reports_dir) && ~isempty(ReportFile)
    bst_report('Export', ReportFile, bst_fullfile(reports_dir, ['report_' ProtocolName '_sources.html']));
end


%% ===== TIME-FREQUENCY =====
% Start a new report (one report for the time-frequency of all the subjects)
bst_report('Start');
% List of conditions to process separately
AllConditions = {'Famous', 'Scrambled', 'Unfamiliar'};
% Channels to display in the screen capture, by order of preference (if the first channel is bad, use the following)
SelChannel = {'EEG070','EEG060','EEG065','EEG050','EEG003'};
% Compute one separate time-frequency average for each subject/run/condition
for iSubj = 1:length(SubjectNames)
    for iRun = 1:6
        % Process: Select data files in: Subject/Run
        sTrialsAll = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname',   SubjectNames{iSubj}, ...
            'condition',     sprintf('run_%02d_sss_notch', iRun));
        % Loop on the conditions
        for iCond = 1:length(AllConditions)
            % Comment describing this average
            strComment = [SubjectNames{iSubj}, ' / ', sprintf('run_%02d', iRun), ' / ', AllConditions{iCond}];
            disp(['BST> ' strComment]);
            % Find the first good channel in the display list
            if isempty(BadChannels{iSubj}{iRun})
                iSel = 1;
            else 
                iSel = find(~ismember(SelChannel,BadChannels{iSubj}{iRun}), 1);
            end
            % Process: Select file comments with tag: Avg
            sTrialsCond = bst_process('CallProcess', 'process_select_tag', sTrialsAll, [], ...
                'tag',    [AllConditions{iCond}, '_trial'], ...
                'search', 1, ...  % Search the file names
                'select', 1);  % Select only the files with the tag
            % Process: Time-frequency (Morlet wavelets), averaged across trials
            sTimefreq = bst_process('CallProcess', 'process_timefreq', sTrialsCond, [], ...
                'sensortypes', 'MEG MAG, EEG', ...
                'edit',        struct(...
                     'Comment',         ['Avg: ' AllConditions{iCond} ', Power, 6-60Hz'], ...
                     'TimeBands',       [], ...
                     'Freqs',           [6, 6.8, 7.6, 8.6, 9.7, 11, 12.4, 14, 15.8, 17.9, 20.2, 22.8, 25.7, 29, 32.7, 37, 41.7, 47.1, 53.2, 60], ...
                     'MorletFc',        1, ...
                     'MorletFwhmTc',    3, ...
                     'ClusterFuncTime', 'none', ...
                     'Measure',         'power', ...
                     'Output',          'average', ...
                     'RemoveEvoked',    0, ...
                     'SaveKernel',      0), ...
                'normalize',   'none');  % None: Save non-standardized time-frequency maps
            % Process: Extract time: [-200ms,900ms]
            sTimefreq = bst_process('CallProcess', 'process_extract_time', sTimefreq, [], ...
                'timewindow', [-0.2, 0.9], ...
                'overwrite',  1);
            % Screen capture of one sensor
            hFigTf = view_timefreq(sTimefreq.FileName, 'SingleSensor', SelChannel{iSel});
            bst_report('Snapshot', hFigTf, strComment, 'Time-frequency', [200, 200, 400, 250]);
            close(hFigTf);
        end
    end
end
% Save report
ReportFile = bst_report('Save', []);
if ~isempty(reports_dir) && ~isempty(ReportFile)
    bst_report('Export', ReportFile, bst_fullfile(reports_dir, ['report_' ProtocolName '_timefreq.html']));
end



