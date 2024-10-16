function fs_segmentation(subjectsDir, sourceSubject)

log_file = fullfile(subjectsDir, sourceSubject, 'fs_segmentation_log.txt');
diary(log_file);

cmd = [...
    'gtmseg --s ' sourceSubject 
    ...];
fs_execute(cmd)

full_cmd = sprintf('export FREESURFER_HOME=%s; source $FREESURFER_HOME/SetUpFreeSurfer.sh; %s', ...
    getenv('FREESURFER_HOME'), cmd);

system(full_cmd);

disp(['Segmentation by gtmseg completed for ', sourceSubject])

diary off
end