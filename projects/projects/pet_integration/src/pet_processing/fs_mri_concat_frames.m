function outputVolume = fs_mri_concat_frames (subjectsDir, sourceSubject, templatePETDirectory, inputVolumes)
cd(subjectsDir);

mriConcatArgs = strjoin(inputVolumes, ' ');
outputVolume = fullfile(templatePETDirectory, 'template_mean.mgz');
 
%--mask;

cmd = ['mri_concat '...
    mriConcatArgs ' --mean ' ...
    '--o ' outputVolume ...
    ];


 fs_execute(cmd)
 outputVolume = 'template_mean.mgz';
end