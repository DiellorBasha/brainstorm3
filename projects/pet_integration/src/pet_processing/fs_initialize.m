
% fs_initialize - Initializes the MATLAB environment for FreeSurfer processing
% This function retrieves environment variables from the system, constructs
% subject directories if they do not exist, and sets MATLAB environment
% variables for FreeSurfer processing.
%
% Inputs:
%   subjectList - A table, structure, or cell array of subject IDs

% Define environment variables within the MATLAB session environment. These
% paths should match exactly those related to your FreeSurfer environment.
% 
subjectsDir = '/export01/data/toolboxes/freesurfer/subjects/';
inputMRIDir = '/export02/data/dbasha/PAD7/data_release_7.0/mri/wave1/';
inputPETDir = '/export02/data/dbasha/PAD7/data_release_7.0/pet/';
inputDerivativesDir = '/export02/data/dbasha/PAD7/data_release_7.0/pet/derivatives/vlpp/';
freesurferHome = '/export01/data/toolboxes/freesurfer/';
stormNetUtilsDir = '/export01/data/dbasha/code/StormNet1/utils';

% Define base output directories (not subject-specific)
outputPETDir = 'pet';
outputPipelineDir = fullfile(outputPETDir, pipelineName);
outputMRIDir = 'anat';
outputSURFDir = 'surf';

% Set MATLAB environment variables for base directories
setenv('SUBJECTS_DIR', subjectsDir);
setenv('INPUT_MRI_DIR', inputMRIDir);
setenv('INPUT_PET_DIR', inputPETDir);
setenv('INPUT_DERIVATIVES_DIR', inputDerivativesDir);
setenv('FREESURFER_HOME', freesurferHome);
setenv('OUTPUT_PIPELINE_DIR', outputPipelineDir);
setenv('STORMNET_UTILS_DIR', stormNetUtilsDir);

% Display the environment variables to verify they are correctly set
fprintf('Environment setup complete:\n');
fprintf('=====================================\n');
fprintf('FREESURFER_HOME: %s\n', freesurferHome);
fprintf('SUBJECTS_DIR: %s\n', subjectsDir);
fprintf('=====================================\n');
fprintf('INPUT_MRI_DIR: %s\n', inputMRIDir);
fprintf('INPUT_PET_DIR: %s\n', inputPETDir);
fprintf('OUTPUT_PIPELINE_DIR: %s\n', outputPipelineDir);
fprintf('=====================================\n');
fprintf('STORMNET_UTILS_DIR: %s\n', stormNetUtilsDir);




