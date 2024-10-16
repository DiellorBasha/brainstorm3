function cmd = fs_recon_all(subject, t1, varargin)
    % fs_recon_all - Run Freesurfer recon-all with optional arguments
    %
    % Usage:
    %   fs_recon_all(subject, t1, 'Name', Value, ...)
    %
    % Inputs:
    %   subject      - Subject ID for Freesurfer
    %   t1           - Path(s) to the T1-weighted image(s), can be a single path or a cell array of paths
    %
    % Optional Name-Value Pairs:
    %   'FLAIR'       - Path to the FLAIR image (default: '')
    %   'T2'          - Path to the T2 image (default: '')
    %   'threads'     - Number of threads to use (default: 1)
    %   'autorecon'   - Stage of recon-all to run (default: '-all')
    %   'subjectDir'  - Directory where the subject folder will be created (default: '')
    %
    % Example:
    %   fs_recon_all('sub001', 'path/to/t1.nii', 'FLAIR', 'path/to/FLAIR.nii', 'threads', 15, 'autorecon', 'autorecon1', 'subjectDir', '/path/to/output')
    %   fs_recon_all('sub001', {'path/to/t1_1.nii', 'path/to/t1_2.nii'}, 'T2', 'path/to/T2.nii', 'threads', 15, 'autorecon', '-autorecon1', 'subjectDir', '/path/to/output')
    
    % Parse required and optional inputs
    (x) || isstring(x));
    addRequired(p, 't1', @(x) ischar(x) || isstring(x) || iscell(x));
    addParameter(p, 'FLAIR', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'T2', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'threads', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'autorecon', '-all', @(x) ischar(x) || isstring(x));
    addParameter(p, 'subjectDir', '', @(x) ischar(x) || isstring(x));
    
    parse(p, subject, t1, varargin{:});
    
    subject = p.Results.subject;
    t1 = p.Results.t1;
    FLAIR = p.Results.FLAIR;
    T2 = p.Results.T2;
    threads = p.Results.threads;
    autorecon = p.Results.autorecon;
    subjectDir = p.Results.subjectDir;
    
    % Construct the t1 input part of the command
    if iscell(t1)
        t1Cmd = strjoin(cellfun(@(x) ['-i ' x], t1, 'UniformOutput', false));
    else
        t1Cmd = ['-i ' t1];
    end
    
    % Determine the pial refinement command
    if ~isempty(T2)
        pialCmd = ['-T2 ' T2 ' -T2pial'];
    elseif ~isempty(FLAIR)
        pialCmd = ['-FLAIR ' FLAIR ' -FLAIRpial'];
    else
        pialCmd = '';
    end
    
    % Construct the recon-all command
    if isempty(subjectDir)
        cmd = sprintf('recon-all -s %s %s %s %s -threads %d', ...
            subject, t1Cmd, pialCmd, autorecon, threads);
    else
        cmd = sprintf('recon-all -s %s %s %s -%s -sd %s -threads %d', ...
            subject, t1Cmd, pialCmd, autorecon, subjectDir, threads);
    end
    
    % Execute the command
   fs_execute(cmd);
    
    % Display completion message
    disp(['Processing completed for ', subject]);
    disp(cmd);cvv b=-