
function registrationFile = fs_pet_coreg (subjectsDir, sourceSubject, templatePETDirectory, templateFile)

% give full path to volume
templatePET = fullfile(templatePETDirectory, templateFile);
    [~, a , ~] =fileparts(templateFile);
    [~, b , ~] =fileparts(a);

registeredPET = fullfile(templatePETDirectory, strcat(b,'.reg.lta'));

cmd = ['mri_coreg'...
    ' --s ' sourceSubject ...
    ' --mov ' templatePET ...
    ' --reg ' registeredPET ...
    ' --threads 15'];

disp([' ' ...
    'Processing completed for ', sourceSubject])


log_file = fullfile(templatePETDirectory, strcat(b, '_fs_pet_coreg_log.txt'));
diary(log_file);

fs_execute(cmd);
  
 diary off 

registrationFile = registeredPET;

end
