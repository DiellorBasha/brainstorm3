function [MriFileGtm, errMsg, regTable] = pet_gtm(PetFile, fwhm, gtmOpts)
% PET_GTM: Geometric Transfer Matrix (Rousset) regional partial volume correction.
%
% Native Brainstorm implementation of the Rousset 1998 GTM: corrects the regional
% mean activity of EVERY region (cortical + subcortical + WM + CSF + cerebellum) for
% PSF spill-over, by inverting the region-to-region spill-over matrix. Unlike
% Muller-Gartner (pet_pvc), it corrects all tissue classes, not just GM. Operates
% directly on the Brainstorm volume grid (no NIfTI round-trip, no SPM segmentation):
% the parcellation comes from the subject's Desikan-Killiany (aparc+aseg) volume.
%
%   W(i,j) = mean over region i of ( PSF (x) 1_region_j )      % spill-over matrix
%   m(i)   = mean of PET over region i                          % observed regional means
%   t      = W \ m                                              % true (corrected) means
%
% The output is a piecewise-constant volume (each region's voxels set to its corrected
% mean), saved as a Brainstorm anatomy node and voxel-aligned with the input PET, so it
% flows through SUVR rescaling and surface projection unchanged.
%
% USAGE:  [MriFileGtm, errMsg, regTable] = pet_gtm(PetFile, fwhm, gtmOpts)
%
% INPUTS:
%   PetFile : static (3D) PET volume in the Brainstorm DB.
%   fwhm    : scanner PSF FWHM in mm (scalar). [] -> from the PET metadata (PET.PsfFwhm, stored
%             at import). Smoothing recorded in PET.SmoothFwhm is added in quadrature
%             (pet_psf_fwhm).
%   gtmOpts : (optional) .AtlasComment (default 'Desikan-Killiany'), .minVox (default 50):
%             labels smaller than minVox voxels are merged into one "brain rest" region.
%
% NOTES:
%   - Voxels where the PET is not finite (outside the PET field of view) are excluded from
%     the regional means and are NaN in the output.
%   - Computing W takes one 3D Gaussian convolution per region: about one minute for the
%     ~100 regions of Desikan-Killiany on a 256^3 grid; finer atlases scale linearly.
%
% OUTPUTS:
%   MriFileGtm : relative path to the corrected (piecewise-constant) volume node.
%   errMsg     : error message, if any.
%   regTable   : struct with .id .name .nvox .observed .corrected (per region). Special ids:
%                -1 extracerebral (or rest, without scalp), -2 air, -3 brain rest (small labels).
%
% SEE ALSO: pet_pvc, pet_scanner_fwhm
%
% Reference: Rousset OG, Ma Y, Evans AC. Correction for partial volume effects in PET:
%            principle and validation. J Nucl Med 1998;39:904-911.

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
% Authors: Diellor Basha, 2026

    MriFileGtm = ''; errMsg = ''; regTable = struct('id',{},'name',{},'nvox',{},'observed',{},'corrected',{});
    if (nargin < 2), fwhm = []; end
    if (nargin < 3) || isempty(gtmOpts), gtmOpts = struct(); end
    if ~isfield(gtmOpts,'AtlasComment') || isempty(gtmOpts.AtlasComment), gtmOpts.AtlasComment = 'Desikan-Killiany'; end
    if ~isfield(gtmOpts,'minVox')       || isempty(gtmOpts.minVox),       gtmOpts.minVox = 50; end

    isProgress = bst_progress('isVisible');
    if ~isProgress, bst_progress('start', 'PET GTM', 'Loading volumes...'); end
    try
        % ----- load PET + parcellation -----
        [sSubject, iSubject] = bst_get('MriFile', PetFile);
        if isempty(sSubject), error('Subject not found for PET file.'); end
        sMriPet = in_mri_bst(PetFile);
        pet = double(sMriPet.Cube(:,:,:,1));
        cubeSize = size(pet);
        voxsize = sMriPet.Voxsize(1:3);
        iAtl = find(strcmp({sSubject.Anatomy.Comment}, gtmOpts.AtlasComment), 1);
        if isempty(iAtl), error('Parcellation "%s" not found.', gtmOpts.AtlasComment); end
        sAtl = in_mri_bst(sSubject.Anatomy(iAtl).FileName);
        L = double(sAtl.Cube(:,:,:,1));
        if ~isequal(size(L), cubeSize), error('Parcellation grid does not match PET grid.'); end

        % ----- PSF -----
        if ~isempty(fwhm) && (any(~isfinite(fwhm(:))) || any(fwhm(:) <= 0))
            error('PSF FWHM must be a positive number (mm).');
        end
        PET = [];
        if isfield(sMriPet, 'PET'), PET = sMriPet.PET; end
        [fwhm, fwhmSrc] = pet_psf_fwhm(PET, fwhm);
        fprintf('BST> PET GTM: PSF FWHM = %.2f mm [%s]\n', fwhm, fwhmSrc);
        fwhm = [fwhm fwhm fwhm];
        isFov = isfinite(pet);                       % outside the PET field of view: excluded

        % ----- regions (complete partition: labels >= minVox, everything else -> "rest" id 0) -----
        bst_progress('text', 'Building GTM regions...');
        ids = unique(L(:)); ids = ids(ids ~= 0);
        keep = ids(arrayfun(@(id) nnz(L==id) >= gtmOpts.minVox, ids));
        small = setdiff(ids, keep);
        if ~isempty(small)                           % small labels -> one "brain rest" region
            keep = [keep(:); -3];
            L(ismember(L, small)) = -3;
        end
        nKeep = numel(keep);
        Lr = zeros(cubeSize);                        % relabelled: 1..nKeep = kept brain regions
        for r = 1:nKeep, Lr(L==keep(r)) = r; end
        restMask = (Lr == 0);
        % Split the leftover ("rest") into an EXTRACEREBRAL nuisance region (skull/scalp/
        % meninges - inside the head, outside the brain) and AIR (outside the head), using
        % the subject's scalp surface. Modeling extracerebral signal as its own region lets
        % GTM remove its spill-in into cortex (important for off-target tracers like tau);
        % a single "rest" region would mix near-zero air with nonzero skull/scalp.
        headMask = local_headmask(sSubject, sMriPet, cubeSize);
        if ~isempty(headMask)
            Lr(restMask &  headMask) = nKeep + 1;    % extracerebral tissue
            Lr(restMask & ~headMask) = nKeep + 2;    % air
            regIds = [keep(:); -1; -2];              % -1 = extracerebral, -2 = air
            R = nKeep + 2;
        else
            Lr(restMask) = nKeep + 1;
            regIds = [keep(:); -1];                  % -1 = rest (no scalp surface available)
            R = nKeep + 1;
            fprintf('BST> PET GTM: no scalp surface -> single rest region (extracerebral not modeled).\n');
        end
        idx = cell(R,1); idxAll = cell(R,1); n = zeros(R,1);
        for i = 1:R
            idxAll{i} = find(Lr==i);                 % all voxels: source of spill-over
            idx{i} = idxAll{i}(isFov(idxAll{i}));    % voxels in the field of view: averaged
            n(i) = numel(idx{i});
        end
        % Drop empty regions (e.g. no air voxels when the head fills the field of view):
        % an empty region gives a 0/0 mean and turns the whole solution into NaN.
        isEmptyReg = (n == 0);
        if any(isEmptyReg)
            idx(isEmptyReg) = []; idxAll(isEmptyReg) = []; n(isEmptyReg) = []; regIds(isEmptyReg) = []; R = numel(n);
        end

        % ----- GTM matrix + observed means -----
        bst_progress('text', sprintf('GTM matrix (%d regions)...', R));
        W = zeros(R,R); m = zeros(R,1);
        for i = 1:R, m(i) = sum(pet(idx{i})) / n(i); end
        for j = 1:R
            mj = zeros(cubeSize); mj(idxAll{j}) = 1;
            hrj = local_gauss3(mj, fwhm, voxsize);
            for i = 1:R, W(i,j) = sum(hrj(idx{i})) / n(i); end
        end

        % ----- solve t = W \ m (robust to ill-conditioning) -----
        c = rcond(W);
        if ~isfinite(c) || (c < 1e-12)
            t = pinv(W) * m;
            fprintf('BST> PET GTM: W ill-conditioned (rcond=%.1e) -> pseudo-inverse.\n', c);
        else
            t = W \ m;
        end

        % ----- piecewise-constant corrected volume -----
        Cout = zeros(cubeSize, 'single');
        for i = 1:R, Cout(idx{i}) = t(i); end
        Cout(~isFov) = NaN;
        sGtm = sMriPet;                              % inherit geometry (same grid -> aligned)
        sGtm.Cube = Cout;
        if isfield(sGtm, 'Histogram'), sGtm.Histogram = []; end   % recomputed on display

        % ----- region table -----
        for i = 1:R
            regTable(i) = struct('id', regIds(i), 'name', local_regname(regIds(i), sAtl, headMask), ...
                                 'nvox', n(i), 'observed', m(i), 'corrected', t(i)); %#ok<AGROW>
        end

        % ----- save as anatomy node -----
        sGtm.Comment = sprintf('%s | PVC GTM %.1fmm', sMriPet.Comment, fwhm(1));
        sGtm = bst_history('add', sGtm, 'gtm', sprintf('Rousset GTM PVC, FWHM=%.2fmm [%s], %d regions, cond=%.1e', ...
                           fwhm(1), fwhmSrc, R, 1/max(c,eps)));
        % Reload the subject: the next call (pet_process, SUVR) looks the new file up in the database
        MriFileGtm = db_add(iSubject, sGtm, 1);
        db_save();
    catch ME
        errMsg = ME.message;
    end
    if ~isProgress, bst_progress('stop'); end
end


function name = local_regname(id, sAtl, headMask)
% Region name from the atlas labels (Labels: {value, name, color} rows), or the special regions.
    switch id
        case -1, if isempty(headMask), name = 'Rest (outside the brain)'; else, name = 'Extracerebral'; end
        case -2, name = 'Air';
        case -3, name = 'Brain rest (small labels)';
        otherwise
            name = sprintf('Label %d', id);
            if isfield(sAtl, 'Labels') && iscell(sAtl.Labels) && (size(sAtl.Labels,2) >= 2)
                iLab = find(cellfun(@(v) isequal(double(v), id), sAtl.Labels(:,1)), 1);
                if ~isempty(iLab), name = sAtl.Labels{iLab,2}; end
            end
    end
end

function headMask = local_headmask(sSubject, sMriRef, cubeSize)
% Voxelize the subject's scalp ("head mask") surface into a filled binary head mask on
% the reference grid. Returns [] if no scalp surface is available.
    headMask = [];
    if ~isfield(sSubject,'Surface') || isempty(sSubject.Surface), return; end
    iScalp = find(strcmpi({sSubject.Surface.SurfaceType}, 'Scalp'), 1);
    if isempty(iScalp), return; end
    try
        tess2mri = tess_interp_mri(sSubject.Surface(iScalp).FileName, sMriRef);
        headMask = logical(tess_mrimask(cubeSize, tess2mri));
    catch
        headMask = [];
    end
end

function vol = local_gauss3(vol, fwhm_mm, voxsize_mm)
    sig = fwhm_mm / 2.35482;
    for d = 1:3
        sv = sig(min(d,numel(sig))) / voxsize_mm(d); r = max(1, ceil(3*sv)); x = -r:r;
        k = exp(-(x.^2)/(2*sv^2)); k = k/sum(k);
        sh = ones(1,3); sh(d) = numel(k);
        vol = convn(vol, reshape(k, sh), 'same');
    end
end
