function [MriFilePvc, errMsg, fileTag] = pet_pvc(PetFile, MriFileRef, fwhm, pvcOpts)
% PET_PVC: Partial volume correction for PET volumes using PETPVE12 (Muller-Gartner method).
%
% USAGE:
%   [MriFilePvc, errMsg, fileTag] = pet_pvc(PetFile, MriFileRef, fwhm)
%   [MriFilePvc, errMsg, fileTag] = pet_pvc(PetFile, MriFileRef, fwhm, pvcOpts)
%
% DESCRIPTION:
%   Performs partial volume correction on a PET volume using the Muller-Gartner
%   method implemented in PETPVE12. The function:
%     1. Exports the reference MRI to NIfTI in a temp directory
%     2. Runs SPM Segment on the MRI to produce probabilistic tissue maps (c1, c2, c3)
%     3. Exports the PET volume to NIfTI in the same temp directory
%     4. Calls PETPVE12's geg_PVEcorrection with the tissue maps and PET
%     5. Imports the corrected PET back into Brainstorm
%     6. Cleans up the temp directory
%
% INPUTS:
%   - PetFile    : PET MRI file to correct (Brainstorm relative path)
%   - MriFileRef : Reference MRI file for tissue segmentation (Brainstorm relative path)
%   - fwhm       : PSF FWHM in mm (scalar or [x y z] vector)
%   - pvcOpts    : (optional) Structure with additional PVC options:
%       .gmThresh    : GM threshold for masking (default: 0.5)
%       .csfZeroing  : If 1, assume CSF signal = 0 (default: 1)
%       .wmcsfMethod : WM/CSF constant estimation method (default: 'threshold')
%                      Options: 'threshold', 'erode'
%       .wmcsfThresh : Threshold/erosion value for WM/CSF estimation (default: 0.5)
%
% OUTPUTS:
%   - MriFilePvc : Relative path to the PVC-corrected PET file in Brainstorm DB
%   - errMsg     : Error message, if any
%   - fileTag    : File tag used for output file
%
% REQUIRES:
%   - SPM12 plugin (for tissue segmentation)
%   - PETPVE12 plugin (for Muller-Gartner PVC)

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

% Initialize outputs
MriFilePvc = [];
errMsg = '';
fileTag = '_pvc';

% ===== PARSE INPUTS =====
if nargin < 3 || isempty(fwhm)
    errMsg = 'PSF FWHM is required for partial volume correction.';
    return;
end
% Expand scalar FWHM to 3D vector
if isscalar(fwhm)
    fwhm = [fwhm fwhm fwhm];
end
% Default PVC options
if nargin < 4 || isempty(pvcOpts)
    pvcOpts = struct();
end
if ~isfield(pvcOpts, 'gmThresh')    || isempty(pvcOpts.gmThresh),    pvcOpts.gmThresh    = 0.5;         end
if ~isfield(pvcOpts, 'csfZeroing')  || isempty(pvcOpts.csfZeroing),  pvcOpts.csfZeroing  = 1;           end
if ~isfield(pvcOpts, 'wmcsfMethod') || isempty(pvcOpts.wmcsfMethod), pvcOpts.wmcsfMethod = 'threshold'; end
if ~isfield(pvcOpts, 'wmcsfThresh') || isempty(pvcOpts.wmcsfThresh), pvcOpts.wmcsfThresh = 0.5;         end

% ===== CHECK PLUGINS =====
% SPM12: Install if needed, then load
[isOk, errInstall] = bst_plugin('Install', 'spm12', 1);
if ~isOk
    errMsg = ['SPM12 plugin is required: ' errInstall];
    return;
end
bst_plugin('Load', 'spm12');
% PETPVE12: Install if needed, then load
[isOk, errInstall] = bst_plugin('Install', 'petpve12', 1);
if ~isOk
    errMsg = ['PETPVE12 plugin is required: ' errInstall];
    return;
end
bst_plugin('Load', 'petpve12');

% ===== PROGRESS BAR =====
isProgress = bst_progress('isVisible');
if ~isProgress
    bst_progress('start', 'PET PVC', 'Loading input volumes...');
end

try
    % ===== LOAD INPUTS =====
    sMriPet = in_mri_bst(PetFile);
    sMriRef = in_mri_bst(MriFileRef);

    % ===== CREATE TEMP DIRECTORY =====
    TmpDir = bst_get('BrainstormTmpDir', 0, 'petpvc');

    % ===== EXPORT REFERENCE MRI TO NIFTI =====
    bst_progress('text', 'Running SPM tissue segmentation...');
    MriNiiFile = bst_fullfile(TmpDir, 'pvc_mri.nii');
    out_mri_nii(sMriRef, MriNiiFile);

    % ===== RUN SPM SEGMENTATION =====
    % Get TPM file from SPM
    TpmFile = bst_get('SpmTpmAtlas');
    if isempty(TpmFile)
        errMsg = 'Could not find SPM TPM atlas.';
        file_delete(TmpDir, 1, 1);
        return;
    end
    % Prepare SPM batch for tissue segmentation
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.preproc.channel(1).vols  = {[MriNiiFile ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg  = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];
    % Tissue classes: GM, WM, CSF (native space only, no warped)
    for iTiss = 1:3
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).tpm    = {[TpmFile, ',' num2str(iTiss)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).ngaus  = [1 1 2];
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).ngaus  = [1, 1, 2];
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).warped = [0 0];
    end
    % Fix ngaus per tissue class
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;  % GM
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;  % WM
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;  % CSF
    % Remaining tissue classes (skull, scalp, air) — needed by SPM but not used
    for iTiss = 4:6
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).tpm    = {[TpmFile, ',' num2str(iTiss)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).ngaus  = [3, 4, 2];
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(iTiss).warped = [0 0];
    end
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    % Warp settings
    matlabbatch{1}.spm.spatial.preproc.warp.mrf     = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup  = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg      = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg   = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm     = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp     = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write    = [0 0];  % No deformation fields needed
    matlabbatch{1}.spm.spatial.preproc.warp.vox      = NaN;
    matlabbatch{1}.spm.spatial.preproc.warp.bb       = [NaN NaN NaN; NaN NaN NaN];
    % Suppress warnings and run
    warning('off', 'MATLAB:RandStream:ActivatingLegacyGenerators');
    spm_jobman('initcfg');
    spm_jobman('run', matlabbatch);
    warning('on', 'MATLAB:RandStream:ActivatingLegacyGenerators');

    % Verify tissue maps were created
    GmFile  = bst_fullfile(TmpDir, 'c1pvc_mri.nii');
    WmFile  = bst_fullfile(TmpDir, 'c2pvc_mri.nii');
    CsfFile = bst_fullfile(TmpDir, 'c3pvc_mri.nii');
    if ~file_exist(GmFile) || ~file_exist(WmFile) || ~file_exist(CsfFile)
        errMsg = 'SPM tissue segmentation failed: tissue maps not found.';
        file_delete(TmpDir, 1, 1);
        return;
    end

    % ===== EXPORT PET TO NIFTI =====
    bst_progress('text', 'Running partial volume correction...');
    PetNiiFile = bst_fullfile(TmpDir, 'pvc_pet.nii');
    out_mri_nii(sMriPet, PetNiiFile);

    % ===== BUILD PETPVE12 JOB =====
    job = struct();
    job.PETdata = {PetNiiFile};
    % Probabilistic tissue maps
    job.SegImgs.Tsegs.tiss1 = {GmFile};
    job.SegImgs.Tsegs.tiss2 = {WmFile};
    job.SegImgs.Tsegs.tiss3 = {CsfFile};
    % PVE correction options
    job.PVEopts.fwhm_PSF = fwhm;
    job.PVEopts.gmthresh = pvcOpts.gmThresh;
    job.PVEopts.TissConv = 0;  % Don't save convolved tissue segments
    % CSF signal handling
    if pvcOpts.csfZeroing
        job.PVEopts.CSFsignal.CSFzeroing = 1;
    else
        job.PVEopts.CSFsignal.CSFcalc.tiss3 = {CsfFile};
    end
    % WM/CSF constant activity estimation
    switch lower(pvcOpts.wmcsfMethod)
        case 'threshold'
            job.PVE_Const_opts.type3.wmcsfthresh = pvcOpts.wmcsfThresh;
        case 'erode'
            job.PVE_Const_opts.type4.erothresh = pvcOpts.wmcsfThresh;
        otherwise
            job.PVE_Const_opts.type3.wmcsfthresh = pvcOpts.wmcsfThresh;
    end

    % ===== RUN PETPVE12 =====
    geg_PVEcorrection(job);

    % ===== IMPORT CORRECTED PET =====
    bst_progress('text', 'Importing corrected PET volume...');
    PvcNiiFile = bst_fullfile(TmpDir, 'pvcpvc_pet.nii');
    if ~file_exist(PvcNiiFile)
        errMsg = 'PETPVE12 correction failed: output file not found.';
        file_delete(TmpDir, 1, 1);
        return;
    end
    % Import corrected PET volume
    sMriPvc = in_mri(PvcNiiFile, 'ALL', 0, 0);
    % Transfer metadata from original PET
    sMriPvc.SCS    = sMriPet.SCS;
    sMriPvc.NCS    = sMriPet.NCS;
    sMriPvc.Header = sMriPet.Header;

    % ===== CLEAN UP TEMP DIRECTORY =====
    file_delete(TmpDir, 1, 1);

    % ===== SAVE IN BRAINSTORM DATABASE =====
    bst_progress('text', 'Saving corrected PET volume...');
    [sSubject, iSubject] = bst_get('MriFile', PetFile);
    if isempty(sSubject) || ~isfield(sSubject, 'Anatomy')
        errMsg = 'Could not find subject for the provided PET file.';
        return;
    end
    % Build output file path
    [folder, base, ext] = fileparts(file_fullpath(PetFile));
    lastUnderscore = find(base == '_', 1, 'last');
    if ~isempty(lastUnderscore)
        newBase = [base(1:lastUnderscore-1), fileTag, base(lastUnderscore:end)];
    else
        newBase = [base, fileTag];
    end
    MriFilePvcFull = file_unique(fullfile(folder, [newBase, ext]));
    MriFilePvc = file_short(MriFilePvcFull);
    % Update comment
    sMriPvc.Comment = file_unique([sMriPet.Comment, fileTag], {sSubject.Anatomy.Comment});
    % Add history entry
    sMriPvc = bst_history('add', sMriPvc, 'pvc', ...
        sprintf('Partial volume correction (Muller-Gartner, FWHM=[%g %g %g]mm)', fwhm(1), fwhm(2), fwhm(3)));
    % Save new MRI in Brainstorm format
    sMriPvc = out_mri_bst(sMriPvc, MriFilePvcFull);
    % Register new MRI in subject
    iAnatomy = length(sSubject.Anatomy) + 1;
    sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
    sSubject.Anatomy(iAnatomy).FileName = MriFilePvc;
    sSubject.Anatomy(iAnatomy).Comment  = sMriPvc.Comment;
    % Update subject structure
    bst_set('Subject', iSubject, sSubject);
    % Refresh tree
    panel_protocols('UpdateNode', 'Subject', iSubject);
    panel_protocols('SelectNode', [], 'anatomy', iSubject, iAnatomy);
    % Save database
    db_save();

catch ME
    errMsg = ME.message;
    % Attempt cleanup on error
    if exist('TmpDir', 'var') && ~isempty(TmpDir) && file_exist(TmpDir)
        file_delete(TmpDir, 1, 1);
    end
end

% Stop progress bar if we started it
if ~isProgress
    bst_progress('stop');
end
