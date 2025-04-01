function [sMri, regionsStats] = mri_mask (MriFile, sMri, AtlasName, RegionName, isBinary)
    % MRI_MASK: Create a mask from a volume atlas (eg. aseg.mgz) and a region name (eg. 'Amygdala L')
    % 
    % USAGE:     [sMri, regionsStats] = mri_mask (MriFile, AtlasName, RegionName)
    %
    % INPUT:
    %    - MriFile   : Full path to the volume atlas (eg. '/path/to/aseg.mgz')
    %    - sMri      : Braistorm MRI structure
    %    - AtlasName : Name of the atlas: {'aseg', 'marsatlas'}
    %    - isBinary  : Return a binary mask (1) or the masked input volume (0)
    % 
    % OUTPUT:
    %    - Labels : sMri
    %    - regionsStats : Struct with the following fields:
    %
    % REFERENCES:

    
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
    % Authors: Diellor Basha, 2025

% tissue fraction calculation

% Parse inputs

% Initialize returned values


% Get subject
if isempty (MriFileSrc)
    sSubject = bst_get('Subject');
else
    [sSubject, iSubject] = bst_get('MriFile', MriFileSrc);
end

% Switch by atlas name
if isempty(Labels) && ~isempty(AtlasName)
    switch lower(AtlasName)
        case 'freesurfer'    % FreeSurfer ASEG + Desikan-Killiany (2006) + Destrieux (2010)
            Labels = mri_getlabels_freesurfer();
            % Get atlas and generate mask
            iAtlas = find(cellfun(@(c) contains(c, 'ASEG'), {sSubject.Anatomy.Comment}));
        case 'marsatlas'     % BrainVISA MarsAtlas (Auzias 2006)
            Labels = mri_getlabels_marsatlas();
        case 'svreg'         % BrainSuite SVREG (Brainsuite1, USCBrain)
            Labels = mri_getlabels_svreg();
    end
end


[sAtlas, ~] = bst_memory('LoadMri', (sSubject.Anatomy(iAtlas).FileName));
if isempty(sAtlas)
    error('Failed to load aseg file.');
end

% Find the index where the second column matches the region name
iLabel = find(cellfun(@(x) contains(x, regionName, 'IgnoreCase', true), Labels(:, 2)));
binaryMask = ismember(sAtlas.Cube, iLabel);
% Flatten the 3D volume to 1D for clustering
voxel_values = sMri.Cube(:);

end