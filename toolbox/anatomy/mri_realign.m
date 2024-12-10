function  [MriFileMean, MriFileAlign, fileTag, sMriMean] = mri_realign (MriFile, Method)
% MRI_REALIGN: Extract frames from dynamic volumes, realign and compute the mean across frames.
%
% USAGE:  [MriFileMean, MriFileAlign, fileTag, sMriMean] = mri_realign(MriFile, Method)
%                         [sMriMean, sMriAlign, fileTag] = mri_realign(sMri, Method)
%         [MriFileMean, MriFileAlign, fileTag, sMriMean] = mri_realign(MriFile)
%                         [sMriMean, sMriAlign, fileTag] = mri_realign(sMri)
%
% INPUTS:
%    - MriFile : Relative path to the Brainstorm Mri file to realign
%    - Method  : Method used for the realignment of the volume (default is spm_realign): 
%                       -'spm_realign' :        uses the SPM plugin  
%                       -'fs_realign'  :        uses Freesurfer 
%
% OUTPUTS:
%    - sMriRealign      : Dynamic Brainstorm Mri structure with realigned frames
%    - sMriMean         : Static Brainstorm Mri structure with mean frame
%    - MriFileMean      : Relative path to the Brainstorm MRI file containing the computed frame mean - static volume (dim4=1) 
%    - MriFileRealign   : Relative path Brainstorm MRI file containing realigned frames - dynamic volume (dim4>1)
%    - errMsg           : Error messages if any
%    - fileTag          : Tag added to the comment/filename

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
% WARRANTY, EXPRESS OR IMPLIED, IN>CLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2016-2023
%          Chinmay Chinara, 2023
%          Diellor Basha, 2024

% ===== LOAD INPUTS =====
% Parse inputs
if (nargin < 2) || isempty(Method)
    Method = 'spm_align';
end
% Progress bar
isProgress = bst_progress('isVisible');
if ~isProgress
    bst_progress('start', 'Realignment', 'Loading input volumes...');
end
    if isstruct(MriFile) % USAGE: [sMriMean, sMriAlign, fileTag] = mri_realign(sMri, Method)
        sMri = MriFile; 
        MriFile = [];
    elseif ischar(MriFile) % USAGE: [MriFileMean, MriFileAlign, fileTag, sMriMean] = mri_realign(MriFile, Method)
        % Get volume in bst format
        sMri = in_mri_bst(MriFile);
    else 
        bst_progress('stop');
        error('Invalid call.');
    end
% Initialize returned variables
    sMriAlign = sMri;
    sMriAlign.Cube=zeros(size(sMri.Cube));
    sMriMean  = [];
    fileTag   = '';
% Define temporary directory for exporting nifti files
TmpDir = bst_get('BrainstormTmpDir', 0, 'mri_frames');
% Initialize output file names
sMriOutNii = bst_fullfile(TmpDir, 'orig.nii'); 
MriFileMean = bst_fullfile(TmpDir, 'meanorig.nii'); % SPM output: static volume with mean of realigned frames
MriFileRealign = bst_fullfile(TmpDir, 'rorig.nii'); % SPM output: dynamic volume with realigned frames
TransfMatFile =  bst_fullfile(TmpDir, 'orig.mat');  % SPM output: transformation matrices for realigned frames

% ====== ALIGN FRAMES =======
numFrames = size(sMri.Cube, 4);  % Number of frames
    if numFrames==1 % If numFrames is 1, volume is static 
        return
    else
        % Remove NaN
        if any(isnan(sMri.Cube(:)))
            sMri.Cube(isnan(sMri.Cube)) = 0;
        end
        out_mri_nii (sMri, sMriOutNii);% Export as Nifti to TmpDir
    end 

switch lower(Method)

    % ===== METHOD: SPM ALIGN =====
    case 'spm_align'
        % Initialize SPM
        [isInstalled, errMsg] = bst_plugin('Install', 'spm12');
        if ~isInstalled
            if ~isProgress
                bst_progress('stop');
            end
            return;
        end
        bst_plugin('SetProgressLogo', 'spm12');

 % === CALL SPM REALIGN ===
  bst_progress('text', sprintf('Aligning %d  frames using SPM Realign...', numFrames));
      % Create realign batch 
    matlabbatch = {};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {TmpDir};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = 'orig';
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
    matlabbatch{2}.spm.util.exp_frames.files(1) = cfg_dep('File Selector (Batch Mode): Selected Files (orig)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{2}.spm.util.exp_frames.frames = Inf;
    matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

        spm('defaults', 'PET');
        spm_jobman('run', matlabbatch);

       % Import the calculated transformation matrices; 
            sTransfMat = in_bst_matrix(TransfMatrix);  
            sMriAligned = in_mri(MriFileRealign, 'ALL', 0, 1);  % Import the realigned dynamic volume            
``          
         % Output file tag
        fileTag = '_spm_align';


    case 'freesurfer'
end

% ===== UPDATE HISTORY ========

sMriAlign.Comment = [sMriAlign.Comment, fileTag]; % Add file tag
    % Add history entry
    sMriAlign = bst_history('add', sMriAlign, 'realign', ['PET Frames realigned using (' Method '): ']);
    
sMriMean.Comment = [fileTag, '_mean']; % Add file tag
    sMriMean = bst_history('add', sMriMean, 'realign', ['PET Frames realigned using (' Method '): ']);
    sMriMean = bst_history('add', sMriMean, 'mean realigned', ['Mean of realigned PET computed using (' Method '): ']);

    file_delete(TmpDir, 1, 1);

% Close progress bar
if ~isProgress
    bst_progress('stop');
end
end