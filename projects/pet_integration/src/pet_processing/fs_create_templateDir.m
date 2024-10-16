function [outputDir, templateDir]= fs_create_templateDir(sourceSubject, inputPET)
%cd (fullfile(subjectsDir, sourceSubject));
% get the current working directory

petDir = fullfile (getenv('SUBJECTS_DIR'), sourceSubject, getenv('OUTPUT_PET_DIR'));

%check if pet folder already exists
if ~exist(petDir, 'dir')
  %Create the pet directory
    mkdir(petDir)
    disp(['Created directory', petDir])
else
    disp(['Directory already exists', petDir])
end

[~,b,~] = fileparts(inputPET);
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
    disp(['Created directory ', outputDir])
else
    disp(['Directory already exists ', outputDir])
end

end