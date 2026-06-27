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
%   fwhm    : PSF FWHM in mm (scalar). [] -> derived from PET metadata (pet_scanner_fwhm).
%   gtmOpts : (optional) .AtlasComment (default 'Desikan-Killiany'), .minVox (default 50).
%
% OUTPUTS:
%   MriFileGtm : relative path to the corrected (piecewise-constant) volume node.
%   errMsg     : error message, if any.
%   regTable   : struct with .id .nvox .observed .corrected (per region).
%
% SEE ALSO: pet_pvc, pet_scanner_fwhm
%
% Reference: Rousset OG, Ma Y, Evans AC. Correction for partial volume effects in PET:
%            principle and validation. J Nucl Med 1998;39:904-911.
%
% Author: Diellor Basha, 2026

    MriFileGtm = ''; errMsg = ''; regTable = struct('id',{},'nvox',{},'observed',{},'corrected',{});
    if (nargin < 2), fwhm = []; end
    if (nargin < 3) || isempty(gtmOpts), gtmOpts = struct(); end
    if ~isfield(gtmOpts,'AtlasComment') || isempty(gtmOpts.AtlasComment), gtmOpts.AtlasComment = 'Desikan-Killiany'; end
    if ~isfield(gtmOpts,'minVox')       || isempty(gtmOpts.minVox),       gtmOpts.minVox = 50; end

    isProgress = bst_progress('isVisible');
    if ~isProgress, bst_progress('start', 'PET GTM', 'Loading volumes...'); end
    try
        % ----- load PET + parcellation -----
        [sSubject, iSubject] = bst_get('MriFile', PetFile);
        if isempty(sSubject), errMsg = 'Subject not found for PET file.'; return; end
        sMriPet = in_mri_bst(PetFile);
        pet = double(sMriPet.Cube(:,:,:,1));
        cubeSize = size(pet);
        voxsize = sMriPet.Voxsize(1:3);
        iAtl = find(strcmp({sSubject.Anatomy.Comment}, gtmOpts.AtlasComment), 1);
        if isempty(iAtl), errMsg = sprintf('Parcellation "%s" not found.', gtmOpts.AtlasComment); return; end
        sAtl = in_mri_bst(sSubject.Anatomy(iAtl).FileName);
        L = double(sAtl.Cube(:,:,:,1));
        if ~isequal(size(L), cubeSize), errMsg = 'Parcellation grid does not match PET grid.'; return; end

        % ----- PSF -----
        if isempty(fwhm)
            PET = []; try w = load(file_fullpath(PetFile),'PET'); if isfield(w,'PET'), PET = w.PET; end; catch; end
            [fwhm, fwhmSrc] = pet_scanner_fwhm(PET);
            fprintf('BST> PET GTM: PSF FWHM = %.1f mm [%s]\n', fwhm, fwhmSrc);
        end
        if isscalar(fwhm), fwhm = [fwhm fwhm fwhm]; end

        % ----- regions (complete partition: labels >= minVox, everything else -> "rest" id 0) -----
        bst_progress('text', 'Building GTM regions...');
        ids = unique(L(:)); ids = ids(ids ~= 0);
        keep = ids(arrayfun(@(id) nnz(L==id) >= gtmOpts.minVox, ids));
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
        idx = cell(R,1); n = zeros(R,1);
        for i = 1:R, idx{i} = find(Lr==i); n(i) = numel(idx{i}); end

        % ----- GTM matrix + observed means -----
        bst_progress('text', sprintf('GTM matrix (%d regions)...', R));
        W = zeros(R,R); m = zeros(R,1);
        for i = 1:R, m(i) = sum(pet(idx{i})) / n(i); end
        for j = 1:R
            hrj = local_gauss3(reshape(double(Lr==j), cubeSize), fwhm, voxsize);
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
        sGtm = sMriPet;                              % inherit geometry (same grid -> aligned)
        sGtm.Cube = Cout;

        % ----- region table -----
        for i = 1:R
            regTable(i) = struct('id', regIds(i), 'nvox', n(i), 'observed', m(i), 'corrected', t(i)); %#ok<AGROW>
        end

        % ----- save as anatomy node -----
        fileTag = '_gtmpvc';
        sGtm.Comment = file_unique([sMriPet.Comment fileTag], {sSubject.Anatomy.Comment});
        sGtm = bst_history('add', sGtm, 'gtm', sprintf('Rousset GTM PVC, FWHM=%gmm, %d regions, cond=%.1e', fwhm(1), R, 1/max(c,eps)));
        [folder, base, ext] = fileparts(file_fullpath(PetFile));
        u = find(base=='_',1,'last');
        if ~isempty(u), newBase = [base(1:u-1) fileTag base(u:end)]; else, newBase = [base fileTag]; end
        MriFileGtmFull = file_unique(fullfile(folder, [newBase ext]));
        MriFileGtm = file_short(MriFileGtmFull);
        out_mri_bst(sGtm, MriFileGtmFull);
        iAnatomy = length(sSubject.Anatomy) + 1;
        sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
        sSubject.Anatomy(iAnatomy).FileName = MriFileGtm;
        sSubject.Anatomy(iAnatomy).Comment  = file_unique([sMriPet.Comment fileTag], {sSubject.Anatomy.Comment});
        bst_set('Subject', iSubject, sSubject);
        panel_protocols('UpdateNode', 'Subject', iSubject);
        db_save();
    catch ME
        errMsg = ME.message;
    end
    if ~isProgress, bst_progress('stop'); end
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
