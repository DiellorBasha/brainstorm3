function fs_T2_pial_refinement(sourceSubject)
% test command that works 
% recon-all -s sub-MTL0018a -T2 $SUBJECTS_DIR/sub-MTL0018/mri/orig/sub-MTL0018_ses-FU48A_T2star.nii.gz -T2pial -autorecon3 -threads 15

% need path to T2 volume
% raw data is in $SUBJECTS_DIR/sourceSubject
% processed data is $SUBJECTS_DIR/targetSubject

subjectsDir = '/export02/export01/data/toolboxes/freesurfer/subjects';

targetSubject = sprintf([sourceSubject 'a']);
log_file = fullfile(subjectsDir, targetSubject, 'recon_T2_refinement_log.txt');

diary(log_file);

anat_folder = fullfile (subjectsDir, sourceSubject, 'mri', 'orig'); % Path to anat folder
    
% List all files in the anat folder
anat_files = dir(fullfile(anat_folder, '*.nii.gz'));
    T2_files = {};
    for j = 1:length(anat_files)
        if contains(anat_files(j).name, '_T2star')
            T2_files{end+1} = anat_files(j).name;
        end
    end

    T2file = fullfile(anat_folder,T2_files{1});

 % before running -T2 commnand, check that it was not already running but
 % crashed
 % freesurfer likes $SUBJECTs_DIR - in system(cmd), use $SUBJECTS_DIR

cmd = ['recon-all -s ' targetSubject ' -T2 ' T2file ' -T2pial -autorecon3 -threads 15'];

full_cmd = sprintf('export FREESURFER_HOME=%s; source $FREESURFER_HOME/SetUpFreeSurfer.sh; %s', ...
                            getenv('FREESURFER_HOME'), cmd);

system(full_cmd);

% if status ~= 0
%     error('Error in refining the cortical surface: %s', cmdout);
% end
     disp(['Processing completed for ', sourceSubject])
     diary off
end
