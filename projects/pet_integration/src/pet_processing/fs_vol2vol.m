function outputVolume = fs_vol2vol (subjectsDir, sourceSubject, templatePETDirectory, templateFile, inputRegistrationFile)

cd(subjectsDir);
templatePET = fullfile(templatePETDirectory, templateFile);

[~, a , ~] =fileparts(templateFile);
[~, b , ~] =fileparts(a);

T1Volume = fullfile(subjectsDir, sourceSubject, '/mri/T1.mgz');

outputVolume = fullfile(templatePETDirectory, strcat(b, '_coreg.mgz'));

cmd = ['mri_vol2vol '...
                      ' --mov ' templatePET ...
                      ' --targ ' T1Volume ...
                      ' --lta ' inputRegistrationFile ...
                      ' --o '   outputVolume ...
                ];

fs_execute(cmd)

outputVolume = outputVolume;
end