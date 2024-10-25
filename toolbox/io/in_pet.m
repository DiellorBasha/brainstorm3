function [PET, vox2ras, tReorient] = in_pet(PetFile, FileFormat, isInteractive, isNormalize)
% IN_pet: Detect file format and load PET file.
% 
% USAGE:  in_pet(PetFile, FileFormat='ALL', isInteractive=1, isNormalize=0)
% INPUT:
%     - PetFile       : full path to a PET file
%     - FileFormat    : Format of the input file (default = 'ALL')
%     - isInteractive : 0 or 1
%     - isNormalize   : If 1, converts values to uint8 and scales between 0 and 1
% OUTPUT:
%     - PET       : Standard brainstorm structure for PET volumes
%     - vox2ras   : [4x4] transformation matrix: voxels 0-based to RAS coordinates
%                   (corresponds to MNI coordinates if the volume is registered to the MNI space)
%     - tReorient : [4x4] transformation matrix: (voxels 0-based scanner) TO (voxels 0-based Brainstorm)

% NOTES:
%     - PET structure:
%         |- Voxsize:   [x y z], size of each PET voxel, in millimeters
%         |- Cube:      PET volume 
%         |- SCS:       Subject Coordinate System definition
%         |    |- NAS:    [x y z] coordinates of the nasion, in voxels
%         |    |- LPA:    [x y z] coordinates of the left pre-auricular point, in voxels
%         |    |- RPA:    [x y z] coordinates of the right pre-auricular point, in voxels
%         |    |- R:      Rotation to convert PET coordinates -> SCS coordinatesd
%         |    |- T:      Translation to convert PET coordinates -> SCS coordinatesd
%         |- Landmarks[]: Additional user-defined landmarks
%         |- NCS:         Normalized Coordinates System (Talairach, MNI, ...)     
%         |    |- AC:             [x y z] coordinates of the Anterior Commissure, in voxels
%         |    |- PC:             [x y z] coordinates of the Posterior Commissure, in voxels
%         |    |- IH:             [x y z] coordinates of any Inter-Hemispheric point, in voxels
%         |- Comment:     PET description

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2008-2023

% Parse inputs
if (nargin < 4) || isempty(isNormalize)
    isNormalize = 0;
end
if (nargin < 3) || isempty(isInteractive)
    isInteractive = 1;
end
if (nargin < 2) || isempty(FileFormat)
    FileFormat = 'ALL';
end
% Get current byte order
ByteOrder = bst_get('ByteOrder');
if isempty(ByteOrder)
    ByteOrder = 'n';
end
% Initialize returned variables
PET = [];
vox2ras = [];
tReorient = [];

% ===== GUNZIP FILE =====
TmpDir = [];
if ~iscell(PetFile)
    % Get file extension
    [fPath, fBase, fExt] = bst_fileparts(PetFile);
    % If file is gzipped
    if strcmpi(fExt, '.gz')
        % Get temporary folder
        TmpDir = bst_get('BrainstormTmpDir', 0, 'importpet');
        % Target file
        gunzippedFile = bst_fullfile(TmpDir, fBase);
        % Unzip file
        res = org.brainstorm.file.Unpack.gunzip(PetFile, gunzippedFile);
        if ~res
            error(['Could not gunzip file "' PetFile '" to:' 10 gunzippedFile ]);
        end
        % Import gunzipped file
        PetFile = gunzippedFile;
        [fPathTmp, fBase, fExt] = bst_fileparts(PetFile);
    end
    % Default comment
    Comment = fBase;
else
    Comment = 'PET';
    fBase = [];
    fPath = [];
end

                
%% ===== DETECT FILE FORMAT =====
isMni = ismember(FileFormat, {'ALL-MNI', 'ALL-MNI-ATLAS'});
isAtlas = ismember(FileFormat, {'ALL-ATLAS', 'ALL-MNI-ATLAS', 'SPM-TPM'});
if ismember(FileFormat, {'ALL', 'ALL-ATLAS', 'ALL-MNI', 'ALL-MNI-ATLAS'})
    % Switch between file extensions
    switch (lower(fExt))
        % case {'.ima', '.dim'},        FileFormat = 'GIS';
        case {'.img','.hdr','.nii'},  FileFormat = 'Nifti1';
        % case {'.mgz','.mgh'},         FileFormat = 'MGH';
        % case {'.mnc','.mni'},         FileFormat = 'MINC';
        case '.mat',                  FileFormat = 'BST';
        otherwise,                    error('File format could not be detected, please specify a file format.');
    end
end

% ===== LOAD PET =====
% Switch between file formats
switch (FileFormat)   
    % case 'GIS'
    %     PET = in_pet_gis(PetFile, ByteOrder);
    case {'Nifti1', 'Analyze'}
        if isInteractive
            [PET, vox2ras, tReorient] = in_pet_nii(PetFile, 1, [], []);
        else
            [PET, vox2ras, tReorient] = in_pet_nii(PetFile, 1, 1, 0);
        end
    case 'MGH'
        if isInteractive
            [PET, vox2ras, tReorient] = in_pet_mgh(PetFile, [], []);
        else
            petDir = bst_fileparts(PetFile);
            isReconAllClinical = ~isempty(file_find(petDir, 'synthSR.mgz', 2));
            if isReconAllClinical
                [PET, vox2ras, tReorient] = in_pet_mgh(PetFile, 0, 1);
            else
                [PET, vox2ras, tReorient] = in_pet_mgh(PetFile, 1, 0);
            end
        end
    case 'MINC'
        error('Not supported yet');
    case 'BST'
        % Check that the filename contains the 'subjectimage' tag
        if ~isempty(strfind(lower(fBase), 'subjectimage'))
            PET = load(PetFile);
        end
    case 'SPM-TPM'
        error('Not supported yet');
        %PET = in_pet_tpm(PetFile);
    otherwise
        error(['Unknown format: ' FileFormat]);
end
% If nothing was loaded
if isempty(PET)
    return
end
% Default comment: File name
if ~isfield(PET, 'Comment') || isempty(PET.Comment)
    PET.Comment = Comment;
end
% Prepare the history of transformations
if ~isfield(PET, 'InitTransf') || isempty(PET.InitTransf)
    PET.InitTransf = cell(0,2);
end
% If a world/scanner transformation was defined: save it
if ~isempty(vox2ras)
    PET.InitTransf(end+1,[1 2]) = {'vox2ras', vox2ras};
end
% If an automatic reorientation of the volume was performed: save it
if ~isempty(tReorient)
    PET.InitTransf(end+1,[1 2]) = {'reorient', tReorient};
end


%% ===== NORMALIZE VALUES =====
% Remove NaN
if any(isnan(PET.Cube(:)))
    PET.Cube(isnan(PET.Cube)) = 0;
end
% Simplify data type
if ~isa(PET.Cube, 'uint8') && ~isAtlas
    % If only int values between 0 and 255: Reduce storage size by forcing to uint8 
    if (max(PET.Cube(:)) <= 255) && (min(PET.Cube(:)) >= 0) && (max(abs(PET.Cube(:) - round(PET.Cube(:)))) < 1e-10)
        PET.Cube = uint8(PET.Cube);
    % Normalize if the cube is not already in uint8 (and if not loading an atlas)
    elseif isNormalize && ~strcmpi(FileFormat, 'ALL-MNI')
        % Convert to double for calculations
        PET.Cube = double(PET.Cube);
        % Start values at zeros
        PET.Cube = PET.Cube - min(PET.Cube(:));
        % Normalize between 0 and 255 and save as uint8
        PET.Cube = uint8(PET.Cube ./ max(PET.Cube(:)) .* 255);
    end
end


%% ===== CONVERT OLD STRUCTURES TO NEW ONES =====
% Apply a coordinates correction
correction = [.5 .5 0];
if isfield(PET, 'SCS') && isfield(PET.SCS, 'FiducialName') && ~isempty(PET.SCS.FiducialName) && isfield(PET.SCS, 'mmCubeFiducial') && ~isempty(PET.SCS.mmCubeFiducial)
    % === NASION ===
    iNas = find(strcmpi(PET.SCS.FiducialName, 'nasion') | strcmpi(PET.SCS.FiducialName, 'NAS'));
    if ~isempty(iNas)
        PET.SCS.NAS = PET.SCS.mmCubeFiducial(:, iNas)' + correction;
    end
    % === LPA ===
    iLpa = find(strcmpi(PET.SCS.FiducialName, 'LeftPreA') | strcmpi(PET.SCS.FiducialName, 'LPA'));
    if ~isempty(iLpa)
        PET.SCS.LPA = PET.SCS.mmCubeFiducial(:, iLpa)' + correction;
    end
    % === RPA ===
    iRpa = find(strcmpi(PET.SCS.FiducialName, 'RightPreA') | strcmpi(PET.SCS.FiducialName, 'RPA'));
    if ~isempty(iRpa)
        PET.SCS.RPA = PET.SCS.mmCubeFiducial(:, iRpa)' + correction;
    end
    % Remove old fields
    PET.SCS = rmfield(PET.SCS, 'mmCubeFiducial');
    PET.SCS = rmfield(PET.SCS, 'FiducialName');
end
if isfield(PET, 'talCS') && isfield(PET.talCS, 'FiducialName') && ~isempty(PET.talCS.FiducialName) && isfield(PET.talCS, 'mmCubeFiducial') && ~isempty(PET.talCS.mmCubeFiducial)
    NCS = db_template('NCS');
    % === AC ===
    iAc = find(strcmpi(PET.talCS.FiducialName, 'AC'));
    if ~isempty(iAc)
        NCS.AC = PET.talCS.mmCubeFiducial(:, iAc)' + correction;
    end
    % === PC ===
    iPc = find(strcmpi(PET.talCS.FiducialName, 'PC'));
    if ~isempty(iPc)
        NCS.PC = PET.talCS.mmCubeFiducial(:, iPc)' + correction;
    end
    % === IH ===
    iIH = find(strcmpi(PET.talCS.FiducialName, 'IH') | strcmpi(PET.talCS.FiducialName, 'IC'));
    if ~isempty(iIH)
        NCS.IH = PET.talCS.mmCubeFiducial(:, iIH)' + correction;
    end
    % Add new field
    PET.NCS = NCS;
    % Remove old fields
    PET = rmfield(PET, 'talCS');
end


%% ===== READ FIDUCIALS FROM BIDS JSON =====
if ~isempty(fPath)
    % Look for adjacent .json file with fiducials definitions (NAS/LPA/RPA)
    jsonFile = bst_fullfile(fPath, [fBase, '.json']);
    % If json file exists
    if file_exist(jsonFile)
        % Load json file: 0-based voxel coordinates
        json = bst_jsondecode(jsonFile);
        [sFid, msg] = process_import_bids('GetFiducials', json, 'voxel');
        if ~isempty(msg)
            disp(['BIDS> ' jsonFile ': ' msg]);
        end
        % If there are fiducials defined in the json file
        if ~isempty(sFid)
            % Apply re-orientation of the volume to the fiducials coordinates
            iTransf = find(strcmpi(PET.InitTransf(:,1), 'reorient'));
            if ~isempty(iTransf)
                tReorient = PET.InitTransf{iTransf(1),2};  % Voxel 0-based transformation, from original to Brainstorm
                fidNames = fieldnames(sFid);
                for f = fidNames(:)'
                    if ~isempty(sFid.(f{1}))
                        sFid.(f{1}) = (tReorient * [sFid.(f{1}), 1]')';
                        sFid.(f{1}) = sFid.(f{1})(1:3);
                    end
                end
            end
            % Convert from (0-based VOXEL) to (1-based voxel) to (PET)
            if ~isempty(sFid.NAS)
                PET.SCS.NAS = (sFid.NAS + 1) .* PET.Voxsize;
            end
            if ~isempty(sFid.LPA)
                PET.SCS.LPA = (sFid.LPA + 1) .* PET.Voxsize;
            end
            if ~isempty(sFid.RPA)
                PET.SCS.RPA = (sFid.RPA + 1) .* PET.Voxsize;
            end
            if ~isempty(sFid.AC)
                PET.NCS.AC = (sFid.AC + 1) .* PET.Voxsize;
            end
            if ~isempty(sFid.PC)
                PET.NCS.PC = (sFid.PC + 1) .* PET.Voxsize;
            end
            if ~isempty(sFid.IH)
                PET.NCS.IH = (sFid.IH + 1) .* PET.Voxsize;
            end
        end
    end
end


%% ===== COMPUTE SCS TRANSFORMATION =====
% If SCS was defined but transformation not computed
if isfield(PET, 'SCS') && all(isfield(PET.SCS, {'NAS','LPA','RPA'})) ...
                       && ~isempty(PET.SCS.NAS) && ~isempty(PET.SCS.RPA) && ~isempty(PET.SCS.LPA) ...
                       && (~isfield(PET.SCS, 'R') || isempty(PET.SCS.R))
    try
        % Compute transformation
        scsTransf = cs_compute(PET, 'scs');
        % If the SCS fiducials stored in the PET file are not valid : ignore them
        if isempty(scsTransf)
            PET = rmfield(PET, 'SCS');
        % Else: use SCS transform
        else
            PET.SCS.R      = scsTransf.R;
            PET.SCS.T      = scsTransf.T;
            PET.SCS.Origin = scsTransf.Origin;
        end
    catch
        bst_error('Impossible to identify the SCS coordinate system with the specified coordinates.', 'PET Viewer', 0);
    end
end

%% ===== SAVE MNI TRANSFORMATION =====
if isMni && ~isempty(vox2ras) && (~isfield(PET, 'NCS') || ~isfield(PET.NCS, 'R') || isempty(PET.NCS.R))
    % 2nd operation: Change reference from (0,0,0) to (1,1,1)
    vox2ras = vox2ras * [1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];
    % 1st operation: Convert from PET(mm) to voxels
    vox2ras = vox2ras * diag(1 ./ [PET.Voxsize, 1]);
    % Copy MNI transformation to output structure
    PET.NCS.R = vox2ras(1:3,1:3);
    PET.NCS.T = vox2ras(1:3,4);
    % Compute default fiducials positions based on MNI coordinates
    PET = pet_set_default_fid(PET);
end


%% ===== DELETE TEMPORARY FILE =====
if ~isempty(TmpDir)
    file_delete(TmpDir, 1, 1);
end

