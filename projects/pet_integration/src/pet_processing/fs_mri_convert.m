function fs_mri_convert (subjectsDir, sourceSubject, inputPETVolume)

cd (fullfile(subjectsDir, sourceSubject));
% get the current working directory
petDir = fullfile (subjectsDir, sourceSubject, 'pet');

%check if pet folder already exists
if ~exist(petDir, 'dir')
  %Create the pet directory
    mkdir pet
    disp(['Created directory', petDir])
else
    disp(['Directory already exists', petDir])
end

[~,b,~] = fileparts(inputPETVolume);
[~,outputDirectory,~] = fileparts(b);
outputDirectory={outputDirectory};
outputDir = extractRelativePaths(outputDirectory, '_trc-');
outputDir = outputDir{1};

templateDir = fullfile(petDir,outputDir);

%check if pet folder already exists
if ~exist(templateDir, 'dir')
  %Create the pet directory
    cd(petDir)
    mkdir(outputDir)
    disp(['Created directory', outputDir])
else
    disp(['Directory already exists', outputDir])
end

templateFile = fullfile (petDir, outputDir, 'template.nii.gz');
% Create PET template with frame 0 for static PET with one frame
cmd = ['mri_convert '...
    inputPETVolume ' --frame 0 '...
    templateFile];

% Create PET template with averaged frames if you more than 1 frame
cmd = ['mri_concat '...
    inputPETVolume '--mean ' ...
    '--o ' templateFile ...
    ];

full_cmd = sprintf('export FREESURFER_HOME=%s; source $FREESURFER_HOME/SetUpFreeSurfer.sh; %s', ...
    getenv('FREESURFER_HOME'), cmd);
system(full_cmd);
% end

disp(['Template made for ', inputPETVolume])

cd(subjectsDir);
end