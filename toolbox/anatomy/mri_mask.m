function [MriFileMask, errMsg, fileTag, binBrainMask] = mri_roimask(MriFileSrc, MriFileRef, Method, SubMethod, isRescale)
% mri_roimask: Skull stripping on 'MriFileSrc' using 'MriFileRef' as reference MRI.
%                 Both volumes must have the same Cube and Voxel size
%
% USAGE:  [MriFileMask, errMsg, fileTag, binBrainMask] = mri_roimask(MriFileSrc, MriFileRef, Method, SubMethod, isRescale)
%            [sMriMask, errMsg, fileTag, binBrainMask] = mri_roimask(sMriSrc,    sMriRef,    Method, SubMethod, isRescale)
%
%
% INPUTS:
%    - MriFileSrc   : MRI structure or MRI file to apply skull stripping on
%    - MriFileRef   : MRI structure or MRI file to find brain masking for skull stripping
%                     If empty, the Default MRI for that Subject with 'MriFileSrc' is used
%    - Method       : If 'ASEG', use 
%                     If 'DKT', 
%                     If 'Desikan-Killiany'
                      If 'Destrieux'
%    - SubMethod    : 
%    - isRescale    :
% OUTPUTS:
%    - MriFileMask  : MRI structure or MRI file after skull stripping
%    - errMsg       : Error message. Empty if no error
%    - fileTag      : Tag added to the comment and filename
%    - binBrainMask : Volumetric binary mask of the skull stripped 'MriFileRef' reference MRI
%
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

% ===== PARSE INPUTS =====
% Parse inputs
if (nargin < 5)
    isRescale = 0;
end
if (nargin < 4)
    SubMethod = [];
end
if (nargin < 3)
    Method = [];
end

% Initialize outputs
MriFileMask  = [];
errMsg       = '';
fileTag      = '';
binBrainMask = [];
atlasDir = bst_fileparts(file_fullpath(sDefSubject.FileName));
% Return if invalid Method
if isempty(Method) || strcmpi(Method, 'Skip')
    return;
end

% Progress bar
isProgress = bst_progress('isVisible');
if ~isProgress
    bst_progress('start', 'MRI skull stripping', 'Loading input volumes...');
end
% USAGE: mri_reslice(sMriSrc, sMriRef)
if isstruct(MriFileSrc)
    sMriSrc = MriFileSrc;
    sMriRef = MriFileRef;
    MriFileSrc = [];
    MriFileRef = [];
% USAGE: mri_reslice(MriFileSrc, MriFileRef)
elseif ischar(MriFileSrc)
    % Get the default MRI for this subject
    if isempty(MriFileRef)
        sSubject = bst_get('MriFile', MriFileSrc);
        MriFileRef = sSubject.Anatomy(sSubject.iAnatomy).FileName;
    end
    % Load MRI volumes
    sMriSrc = in_mri_bst(MriFileSrc);
    sMriRef = in_mri_bst(MriFileRef);
else
    error('Invalid call.');
end

% Check that same size
refSize = size(sMriRef.Cube(:,:,:,1));
srcSize = size(sMriSrc.Cube(:,:,:,1));
if ~all(refSize == srcSize) || ~all(round(sMriRef.Voxsize(1:3) .* 1000) == round(sMriSrc.Voxsize(1:3) .* 1000))
    errMsg = 'Skull stripping cannot be performed if the reference MRI has different size';
    return
end

% === SKULL STRIPPING ===
% Reset any previous logo
bst_plugin('SetProgressLogo', []);
switch lower(Method)
    case 'aseg'
        % Get subject
        if isempty (MriFileSrc)
            sSubject = bst_get('Subject');
        else
            [sSubject, iSubject] = bst_get('MriFile', MriFileSrc);
        end
        % Get atlas and generate mask
        iAtlas = find(cellfun(@(c) contains(c, 'ASEG'), {sSubject.Anatomy.Comment}));
        [sAtlas, ~] = bst_memory('LoadMri', (sSubject.Anatomy(iAtlas).FileName));
        if isempty(sAtlas)
            error('Failed to load aseg file.');
        end
        % Extract subregion masks
        switch lower(SubMethod)
            case 'brainmask'
                iLabel = 1:255; % cerebrum
            case 'cortex'
                iLabel = [3,42];
            case 'gm'
                iLabel = [3, 8, 9, 10, 11, 12, 13, 17, 18, 19, 26, 27, 42, 47, 48, 49, 50, 51, 52, 52, 54, 55, 58, 59];
            case 'wm'
                iLabel = [2, 7, 41, 46];
            case 'ventricles'
                iLabel = [4, 5, 14, 15, 43, 44, 72];
            case 'hippocampus'
                iLabel = [17, 53];
            case 'brainstem'
                iLabel = 16;
            case 'cerebellum'
                iLabel = [7, 8, 46, 47];
            case 'cerebellum-gm'
                iLabel = [8,47];
        end
          binBrainMask = ismember(sAtlas.Cube, iLabel);
          filesDel = '';
          Method = [SubMethod '-' Method]; % for the file tag
    otherwise
        errMsg = ['Invalid skull stripping method: ' Method];
        return
end
% Reset logo
bst_progress('removeimage');

% Apply brain mask
sMriMask = sMriSrc;
sMriMask.Cube(~binBrainMask) = 0;
% File tag
fileTag = sprintf('_masked_%s', lower(Method));

% === RESCALING ===
    if isRescale
    % rescale the source cube to the average of the masked region
    sMriMask.Cube = double(sMriMask.Cube);  % Convert to double for calculations 
    sMriSrc.Cube = double(sMriSrc.Cube);
    % shift: start values at 1
        %sMriSrc.Cube = (sMriSrc.Cube - min(sMriSrc.Cube(:)))+1;
        maskMean = mean(sMriMask.Cube(binBrainMask));
        sMriSrc.Cube = sMriSrc.Cube./maskMean;  
    sMriMask = sMriSrc; % output rescaled volume)
    % File tag
    fileTag = sprintf('_rescaled_%s', lower(SubMethod));
    end

% ===== SAVE NEW FILE =====
% Add file tag
sMriMask.Comment = [sMriSrc.Comment, fileTag];
% Save output
if ~isempty(MriFileSrc)
    bst_progress('text', 'Saving new file...');
    % Get subject
    [sSubject, iSubject] = bst_get('MriFile', MriFileSrc);
    % Update comment
    sMriMask.Comment = file_unique(sMriMask.Comment, {sSubject.Anatomy.Comment});
    if isRescale
        % Add history entry
    sMriMask = bst_history('add', sMriMask, 'resample', ['Rescaled to "' SubMethod '" using: ' Method]);
    else
    % Add history entry
    sMriMask = bst_history('add', sMriMask, 'resample', ['Skull stripping with "' Method '" using on default file: ' MriFileRef]);
    end
    % Save new file
    MriFileMaskFull = file_unique(strrep(file_fullpath(MriFileSrc), '.mat', [fileTag '.mat']));
    MriFileMask = file_short(MriFileMaskFull);
    % Save new MRI in Brainstorm format
    sMriMask = out_mri_bst(sMriMask, MriFileMaskFull);
    % Register new MRI
    iAnatomy = length(sSubject.Anatomy) + 1;
    sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
    sSubject.Anatomy(iAnatomy).FileName = MriFileMask;
    sSubject.Anatomy(iAnatomy).Comment  = sMriMask.Comment;
    % Update subject structure
    bst_set('Subject', iSubject, sSubject);
    % Refresh tree
    panel_protocols('UpdateNode', 'Subject', iSubject);
    panel_protocols('SelectNode', [], 'anatomy', iSubject, iAnatomy);
    % Save database
    db_save();
else
    % Return output structure
    MriFileMask = sMriMask;
end

% Delete the temporary files
file_delete(filesDel, 1, 1);
% Close progress bar
if ~isProgress
    bst_progress('stop');
end