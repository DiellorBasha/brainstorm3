function [MriFileOut, errMsg, SurfaceFileOut] = pet_process(PetFile, AtlasName, roiName, maskROI, applyMask, doProject, pvcOpts)
% PET_PROCESS: Script PET processing pipeline (PVC, SUVR rescale, masking) with minimal redundant saving.
%
% INPUTS:
%   - PetFile   : PET file path
%   - AtlasFile : Name of anatomical Atlas
%   - roiName   : Name of the ROI for SUVR rescale (string, can be empty)
%   - maskROI   : Name of the ROI for masking (string, can be empty)
%   - applyMask : Logical, true to apply mask, false otherwise
%   - doProject : Logical, true to project PET to surface, false otherwise
%   - pvcOpts   : (optional) Structure controlling PVC / smoothing / SUVR (see pet_pvc.m).
%                  If provided, PVC is applied before SUVR rescaling. Fields:
%                  .method     'gtm' (default, pet_gtm) | 'mg' (Muller-Gartner, PETPVE12) |
%                              'none' (SKIP PVC).
%                  .fwhm       scanner PSF FWHM in mm for PVC ([] or omitted: stored at import
%                              from the scanner metadata, 6 mm if unknown). Recorded smoothing
%                              is added in quadrature (pet_psf_fwhm).
%                  .SmoothFWHM Gaussian volume smoothing FWHM (mm) applied AFTER PVC and
%                              BEFORE SUVR (default 0). Recorded in PET.SmoothFwhm.
%                  .SuvrOpts   struct passed to pet_suvr. Default: robust reference (1-voxel
%                              erosion + 10% trimmed mean). struct('Reference','plain') gives
%                              the previous behaviour (plain mean over the region, mri_rescale).
%                  The intermediate PVC volume is kept in the database ("| PVC GTM <fwhm>mm").
%
%                  VLPP-STYLE (non-PVC, smoothed) run: pvcOpts = struct('method','none',
%                  'SmoothFWHM',6, 'SuvrOpts',struct('Reference','plain')). This reproduces
%                  the Villeneuve Lab PET pipeline (smooth + plain-reference SUVR, no PVC):
%                  https://github.com/villeneuvelab/vlpp
%
% OUTPUTS:
%   - MriFileOut    : Output MRI file path (string)
%   - errMsg        : Error message, if any
%   - SurfaceFileOut: Output surface file path (string, if projected)
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
MriFileOut = '';
SurfaceFileOut = '';
errMsg = '';
if nargin < 6 || isempty(doProject)
    doProject = 0;
end
if nargin < 7
    pvcOpts = [];
end

try
    % Get Subject for PET file
    [sSubject, iSubject] = bst_get('MriFile', PetFile);
    % Load Atlas file
    [~, iAtlas] = ismember(AtlasName, {sSubject.Anatomy.Comment});
    if iAtlas
        AtlasName = sSubject.Anatomy(iAtlas).FileName;
    end
    sAtlas = in_mri_bst(AtlasName);
    % Load PET file in sMRI structure
    sMri = in_mri_bst(PetFile);
    orgComment = sMri.Comment;

    % --- Partial Volume Correction (before SUVR) ---
    % pvcOpts.method='none' SKIPS PVC (e.g. a VLPP-style non-PVC smoothed pipeline, see header
    % + https://github.com/villeneuvelab/vlpp). 'gtm' (default) or 'mg': pet_pvc.
    doPvc = ~isempty(pvcOpts) && ~(isfield(pvcOpts,'method') && strcmpi(pvcOpts.method,'none'));
    if doPvc
        % Get reference MRI for tissue segmentation
        if isfield(sSubject, 'iAnatomy') && ~isempty(sSubject.iAnatomy)
            MriFileRef = sSubject.Anatomy(sSubject.iAnatomy).FileName;
        else
            MriFileRef = sSubject.Anatomy(1).FileName;
        end
        % FWHM: explicit if provided, else [] -> pet_pvc auto-derives from scanner/metadata.
        if isfield(pvcOpts, 'fwhm'), pvcFwhm = pvcOpts.fwhm; else, pvcFwhm = []; end
        % Run PVC - saves corrected PET to database, returns new file path
        [PvcFile, errMsgPvc] = pet_pvc(PetFile, MriFileRef, pvcFwhm, pvcOpts);
        if ~isempty(errMsgPvc)
            errMsg = ['PVC failed: ' errMsgPvc];
            return;
        end
        % Continue processing with the PVC-corrected file
        PetFile = PvcFile;
        sMri = in_mri_bst(PetFile);
        orgComment = sMri.Comment;
    end

    % --- Optional Gaussian volume smoothing BEFORE SUVR (VLPP-style; default off) ---
    smoothTag = '';
    if ~isempty(pvcOpts) && isfield(pvcOpts,'SmoothFWHM') && ~isempty(pvcOpts.SmoothFWHM) && any(pvcOpts.SmoothFWHM(:) > 0)
        sMri.Cube = local_gauss3(double(sMri.Cube(:,:,:,1)), pvcOpts.SmoothFWHM, sMri.Voxsize);
        sMri = bst_history('add', sMri, 'smooth', sprintf('Gaussian volume smoothing FWHM=%g mm', pvcOpts.SmoothFWHM(1)));
        % Record the smoothing so that a later PVC of this volume adds it to the PSF
        if isfield(sMri, 'PET') && isstruct(sMri.PET)
            if ~isfield(sMri.PET, 'SmoothFwhm'), sMri.PET.SmoothFwhm = []; end
            sMri.PET.SmoothFwhm(end+1) = pvcOpts.SmoothFWHM(1);
        end
        smoothTag = sprintf('_smooth%g', pvcOpts.SmoothFWHM(1));
    end

    % --- SUVR Rescale (robust reference: erosion + trimmed mean, via pet_suvr) ---
    if ~isempty(roiName)
        % Resolve the reference ROI to a binary mask (same atlas/region path as before),
        % then normalize via pet_suvr. Default = eroded + trimmed mean; pvcOpts.SuvrOpts can
        % override (e.g. Erode=0, Robust='mean' for a plain VLPP-style cerebellar reference).
        [~, ~, errMsgMask, ~, binMask] = mri_mask(sMri, sAtlas, roiName, 1);
        if ~isempty(errMsgMask)
            errMsg = errMsgMask;
            return;
        end
        suvrOpts = struct('RefMask', binMask);
        if ~isempty(pvcOpts) && isfield(pvcOpts,'SuvrOpts') && isstruct(pvcOpts.SuvrOpts)
            sf = fieldnames(pvcOpts.SuvrOpts); for ii=1:numel(sf), suvrOpts.(sf{ii}) = pvcOpts.SuvrOpts.(sf{ii}); end
        end
        [sMri, ~] = pet_suvr(sMri, [], suvrOpts);
        fileTag = [smoothTag '_suvr'];
    else
        fileTag = smoothTag;
    end

    % --- Masking (if requested) ---
    if applyMask && ~isempty(maskROI)
        [~, sMriMasked, errMsgMask, fileTagMask] = mri_mask(sMri, sAtlas, maskROI, 1);
        if ~isempty(errMsgMask)
            errMsg = errMsgMask;
            return;
        end
        sMri = sMriMasked;
        % Combine tags for output file
        fileTag = [fileTag, fileTagMask];
    end

    % --- Save output file manually (like mri_realign) ---
    % Insert fileTag before last underscore
    [folder, base, ext] = fileparts(file_fullpath(PetFile));
    lastUnderscore = find(base == '_', 1, 'last');
    if ~isempty(lastUnderscore)
        newBase = [base(1:lastUnderscore-1), fileTag, base(lastUnderscore:end)];
    else
        newBase = [base, fileTag];
    end
    MriFileOutFull = file_unique(fullfile(folder, [newBase, ext]));
    MriFileOut = file_short(MriFileOutFull);

    % Update comment to be unique
    sSubject = bst_get('Subject', iSubject);
    sMri.Comment = file_unique([orgComment, fileTag], {sSubject.Anatomy.Comment});

    % Add history entry
    if ~isempty(roiName)
        sMri = bst_history('add', sMri, 'rescale', sprintf('Rescaled with "%s" (%s)', sAtlas.Comment, roiName));
    end
    if applyMask && ~isempty(maskROI)
        sMri = bst_history('add', sMri, 'mask', sprintf('Masked with "%s" (%s)', sAtlas.Comment, maskROI));
    end

    % Save new MRI in Brainstorm format
    sMri = out_mri_bst(sMri, MriFileOutFull);

    % Register new MRI in subject
    iAnatomy = length(sSubject.Anatomy) + 1;
    sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
    sSubject.Anatomy(iAnatomy).FileName = MriFileOut;
    sSubject.Anatomy(iAnatomy).Comment  = sMri.Comment;
    bst_set('Subject', iSubject, sSubject);

    % Refresh tree and save database
    panel_protocols('UpdateNode', 'Subject', iSubject);
    panel_protocols('SelectNode', [], 'anatomy', iSubject, iAnatomy);
    db_save();

    % --- Overlay processed PET on subject's default MRI ---
    if isfield(sSubject, 'iAnatomy') && ~isempty(sSubject.iAnatomy) && ...
            isfield(sSubject, 'Anatomy') && length(sSubject.Anatomy) >= sSubject.iAnatomy
        refMriFile = sSubject.Anatomy(sSubject.iAnatomy).FileName;
    else
        refMriFile = sSubject.Anatomy(1).FileName;
    end
    view_mri(refMriFile, MriFileOut);

    % --- Project to surface if requested ---
    if doProject
        % Use the same reference MRI as above
        % Use the subject's name as the condition
        Condition = 'PET';
        DisplayUnits = '';
        ProjFrac = [0.1 0.4 0.5];
        [SurfaceFileOut, errProj] = mri_interp_vol2tess(MriFileOut, refMriFile, Condition, DisplayUnits, ProjFrac);
        if ~isempty(errProj)
            errMsg = ['PET processed, but projection failed: ', errProj];
        end
    end

catch ME
    errMsg = ME.message;
end
end

% ===== 3D separable Gaussian volume smoothing (FWHM mm), toolbox-free (mirrors pet_gtm) =====
function vol = local_gauss3(vol, fwhm_mm, voxsize_mm)
    if isscalar(fwhm_mm),  fwhm_mm  = fwhm_mm  * [1 1 1]; end
    if isempty(voxsize_mm) || any(voxsize_mm == 0), voxsize_mm = [1 1 1]; end
    sig = fwhm_mm / 2.35482;
    for d = 1:3
        sv = sig(min(d,numel(sig))) / voxsize_mm(d); r = max(1, ceil(3*sv)); x = -r:r;
        k = exp(-(x.^2)/(2*sv^2)); k = k/sum(k);
        sh = ones(1,3); sh(d) = numel(k);
        vol = convn(vol, reshape(k, sh), 'same');
    end
end