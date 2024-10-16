function volumeFiles = fs_findInputVolume(sourceSubject, sourceSession, stringInFile, varargin)
    % fs_findInputVolume - Find input volume files matching specified criteria
    %
    % Syntax: volumeFiles = fs_findInputVolume(sourceSubject, sourceSession, stringInFile, varargin)
    %
    % Inputs:
    %    sourceSubject - Subject identifier (string)
    %    sourceSession - Session identifier (string)
    %    stringInFile - String to match within file names (string)
    %    varargin - Additional options (name-value pairs)
    %               'SourceSubmodality' - Submodality identifier (string, optional)
    %               'PipelineName' - Name of the processing pipeline (string, optional)
    %               'BaseDir' - Base directory for input files (string, optional)
    %               'FileExtension' - File extension to look for (string, default: '*.nii.gz')
    %
    % Outputs:
    %    volumeFiles - Struct array with fields 'name' and 'path' for matching files

    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'SourceSubmodality', '', @ischar);
    addParameter(p, 'PipelineName', '', @ischar);
    addParameter(p, 'BaseDir', getenv('INPUT_MRI_DIR'), @ischar);
    addParameter(p, 'FileExtension', '*.nii.gz', @ischar);
    parse(p, varargin{:});

    sourceSubmodality = p.Results.SourceSubmodality;
    pipelineName = p.Results.PipelineName;
    baseDir = p.Results.BaseDir;
    fileExtension = p.Results.FileExtension;

    % Construct submodality directory
    if isempty(sourceSubmodality)
        searchDir = fullfile(baseDir, sourceSubject, sourceSession);
    else
        searchDir = fullfile(baseDir, sourceSubject, sourceSession, sourceSubmodality);
    end

    % Debug: Display the constructed directory path
    fprintf('Search directory: %s\n', searchDir);

    % Check if search directory exists
    if ~exist(searchDir, 'dir')
        error('Search directory does not exist: %s', searchDir);
    end

    % List all files in the search directory
    searchFiles = dir(fullfile(searchDir, fileExtension));

    % Debug: Display the number of files found
    fprintf('Number of files found: %d\n', length(searchFiles));

    % Check if any files are found
    if isempty(searchFiles)
        warning('No files found in directory: %s', searchDir);
    end

    % Initialize output struct
    volumeFiles = struct('name', {}, 'path', {});

    % Match files containing the specified string
    volumeTypeIn = strcat(stringInFile);

    % Debug: Display the string to match in filenames
    fprintf('String to match in filenames: %s\n', volumeTypeIn);

    for j = 1:length(searchFiles)
        if contains(searchFiles(j).name, volumeTypeIn)
            volumeFiles(end+1).name = searchFiles(j).name; %#ok<*AGROW>
            volumeFiles(end).path = fullfile(searchFiles(j).folder, searchFiles(j).name);

            % Debug: Display the matched file information
            fprintf('Matched file: %s\n', searchFiles(j).name);
            fprintf('Matched file path: %s\n', fullfile(searchFiles(j).folder, searchFiles(j).name));
        end
    end
end
