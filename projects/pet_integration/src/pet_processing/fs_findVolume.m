function [filePaths, volumeFiles] = fs_findVolume(searchDir, dataType)
    % List all files in the search directory
    % dataType = 'T1w'; 'T2', 'FLAIR', 'pet', 'meg'
    searchFiles = dir(fullfile(searchDir, strcat('*_', dataType, '.nii.gz')));

    % Debug: Display the number of files found
    fprintf('Number of files found: %d\n', length(searchFiles));

    % Check if any files are found
    if isempty(searchFiles)
        warning('No files found in directory: %s', searchDir);
    end

    % Initialize output struct and cell array
    volumeFiles = struct('name', {}, 'path', {});
    filePaths = {};

    % Debug: Display the string to match in filenames
    fprintf('String to match in filenames: %s\n', dataType);

    for j = 1:length(searchFiles)
        if contains(searchFiles(j).name, dataType)
            volumeFiles(end+1).name = searchFiles(j).name; 
            volumeFiles(end).path = fullfile(searchFiles(j).folder, searchFiles(j).name);
            filePaths{end+1} = volumeFiles(end).path; % Add path to cell array
            
            % Debug: Display the matched file information
            fprintf('Matched file: %s\n', searchFiles(j).name);
            fprintf('Matched file path: %s\n', fullfile(searchFiles(j).folder, searchFiles(j).name));
        end
    end
end
