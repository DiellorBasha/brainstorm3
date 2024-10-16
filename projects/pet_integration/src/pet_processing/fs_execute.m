function fs_execute(cmd)

full_cmd = sprintf('export FREESURFER_HOME=%s; source $FREESURFER_HOME/SetUpFreeSurfer.sh; %s', ...
    getenv('FREESURFER_HOME'), cmd);
system(full_cmd)

end