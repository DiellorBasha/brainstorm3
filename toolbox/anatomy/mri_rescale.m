function [MriFileMask, errMsg, fileTag, binBrainMask] = mri_rescale(MriFileSrc, MriFileRef, Method)

% MRI_SKULLSTRIP: Skull stripping on 'MriFileSrc' using 'MriFileRef' as reference MRI.
%                 Both volumes must have the same Cube and Voxel size
%
% USAGE:  [MriFileMask, errMsg, fileTag, binBrainMask] = mri_rescale(MriFileSrc, MriFileRef, Method, subMethod)
%            [sMriMask, errMsg, fileTag, binBrainMask] = mri_rescale(sMriSrc,    sMriRef,    Method, subMethod)
%
% INPUTS:
%    - MriFileSrc   : MRI structure or MRI file to apply skull stripping on
%    - MriFileRef   : MRI structure or MRI file to find brain masking for skull stripping
%                     If empty, the Default MRI for that Subject with 'MriFileSrc' is used
%    - Method       : Subregion from Freesurfer's Color Lookup table
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
% Authors: Diellor Basha, 2024
%          Raymundo Cassani, 2024
%          Chinmay Chinara, 2024
%

% ===== PARSE INPUTS =====
% Parse inputs
if (nargin < 3)
    Method = [];
end

% Initialize outputs
MriFileMask  = [];
errMsg       = '';
fileTag      = '';
binBrainMask = [];

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

% === RESCALING ===

% Reset any previous logo
bst_plugin('SetProgressLogo', []);

[sMriMask, errMsg, fileTag, binBrainMask] = mri_skullstrip(sMriSrc, sMriRef, 'ASEG', Method);

% take average of the masked region and rescale
maskMean = mean(sMriMask.Cube(binBrainMask));
sMriSrc.Cube = sMriSrc.Cube./maskMean;
