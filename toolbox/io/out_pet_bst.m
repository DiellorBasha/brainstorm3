function PET = out_pet_bst( PET, PetFile, Version)
% OUT_PET_BST: Save a Brainstorm PET structure.
% 
% USAGE:  PET = out_pet_bst( PET, PetFile )
%
% INPUT: 
%     - PET     : Brainstorm PET structure
%     - PetFile : full path to file where to save the PET in brainstorm format
%     - Version : 'v6', fastest option, bigger files, no files >2Gb
%                 'v7', slower option, compressed, no files >2Gb (default)
%                 'v7.3', much slower, compressed, allows files >2Gb
% OUTPUT:
%     - PET : Modified PET structure
%
% NOTES:
%     - PET structure:
%         |- Cube:      PET volume 
%         |- Voxsize:   [x y z], size of each PET voxel, in millimeters'
%         |- Tracer:    The name of the PET tracer
%         |- Frames:    The total number of PET frames in the volume
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
% Authors: Francois Tadel, 2008-2012

if nargin < 3
    Version = 'v7';
end

% ===== Clean-up PET structure =====
% Remove (useless or old fieldnames)
Fields2BDeleted     = {'Origin','sag','ax','cor','hFiducials','header','filename'};
SCSFields2BDeleted  = {'Origin','Comment'};

for k = 1:length(Fields2BDeleted)
    if isfield(PET,Fields2BDeleted{k})
        PET = rmfield(PET,Fields2BDeleted{k});
    end
end

if isfield(PET,'SCS')
    for k = 1:length(SCSFields2BDeleted)
        if isfield(PET.SCS,SCSFields2BDeleted{k})
            PET.SCS = rmfield(PET.SCS,SCSFields2BDeleted{k});
        end
    end
end

if isfield(PET,'Landmarks') && ~isempty(PET.Landmarks)
    tmpLandmarks = PET.Landmarks;
    nLandmarks2Remove = 0;
    if isfield(PET,'SCS')
        nLandmarks2Remove = 4; %Remove SCS fiducials from Landmark list
    end
    if isfield(PET,'talCS')
        nLandmarks2Remove = nLandmarks2Remove + 3; % Remove TAL/MNI fiducials
    end
    if isfield(PET.Landmarks, 'Names') && ~isempty(PET.Landmarks.Names)
        PET.Landmarks.Names    = PET.Landmarks.Names(nLandmarks2Remove+1:end);
    end
    if isfield(PET.Landmarks, 'PETmmXYZ') && ~isempty(PET.Landmarks.PETmmXYZ)
        PET.Landmarks.PETmmXYZ = PET.Landmarks.PETmmXYZ(:,nLandmarks2Remove+1:end);
    end
    if isfield(PET.Landmarks, 'Handles')
        PET.Landmarks = rmfield(PET.Landmarks,'Handles');
    end
end

if isfield(PET,'SCS2Landmarks') % don't need to be stored in file
    tmpSCS2Landmarks = PET.SCS2Landmarks;
    PET = rmfield(PET,'SCS2Landmarks');
end

% Remove FileName field
if isfield(PET,'FileName')
    PET = rmfield(PET, 'FileName');
end


% SAVE .mat file
try
    bst_save(PetFile, PET, Version);
catch

end
end


