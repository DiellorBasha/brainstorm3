function [sMriSuvr, info] = pet_suvr(sMriPet, sAseg, Opts)
% PET_SUVR: SUVR normalization with a robust, erosion-cleaned reference region.
%
% Upgrades the plain "mean over an ROI" rescale (mri_rescale) with the two things a reference
% region needs to be trustworthy:
%   - EROSION: peel Opts.Erode voxels off the reference mask so its boundary (which is partial-
%     volume-contaminated by neighbouring tissue/CSF) does not bias the reference downward.
%   - ROBUST estimator: a trimmed mean (default) over the eroded interior, rejecting the residual
%     outlier voxels (mis-segmentation, vessels, spill) that a plain mean is sensitive to.
% SUVR = PET / reference. Reference defaults to the cerebellar cortex (ASEG labels 8 + 47), the
% standard amyloid reference (Centiloid-compatible). For tau the inferior cerebellar GM (SUIT) is
% preferable; pass Opts.RefLabels for the available segmentation.
%
% Both volumes must be on the same grid (PET resliced to the anatomy; ASEG in anatomy space).
%
% USAGE:  [sMriSuvr, info] = pet_suvr(sMriPet, sAseg, Opts)
%
% INPUTS:
%   sMriPet : (PVC'd) static PET MRI struct, resliced to the anatomy grid.
%   sAseg   : ASEG volume atlas struct (.Cube of integer labels) on the same grid. May be []
%             if Opts.RefMask is supplied.
%   Opts    : .RefMask (precomputed binary reference mask; overrides RefLabels/sAseg)
%             .RefLabels (default [8 47]) .Erode (voxels, 1) .Robust ('trim'|'mean'|'median')
%             .TrimPct (each-tail fraction for 'trim', 0.10).
%
% OUTPUTS:
%   sMriSuvr : SUVR volume (Cube ./ reference), Comment/History updated.
%   info     : struct(.RefValue,.nVoxMask,.nVoxEroded,.RefMean,.RefMedian,.RefTrim,.Erode,.Robust).
%
% SEE ALSO: pet_process, pet_pvc, mri_rescale
%
% Author: Diellor Basha, 2026

    if (nargin<3)||isempty(Opts), Opts=struct(); end
    Def=struct('RefMask',[],'RefLabels',[8 47],'Erode',1,'Robust','trim','TrimPct',0.10);
    fn=fieldnames(Def); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=Def.(fn{i}); end; end

    cube = double(sMriPet.Cube(:,:,:,1));
    % Reference mask: precomputed (any region, from the caller) or built from ASEG labels.
    if ~isempty(Opts.RefMask)
        mask = logical(Opts.RefMask);
    elseif ~isempty(sAseg)
        mask = ismember(sAseg.Cube, Opts.RefLabels);
    else
        error('pet_suvr:ref', 'Provide Opts.RefMask, or sAseg + Opts.RefLabels.');
    end
    if ~isequal(size(cube), size(mask))
        error('pet_suvr:grid', 'PET grid %s != reference mask grid %s (reslice PET to anatomy first).', ...
              mat2str(size(cube)), mat2str(size(mask)));
    end
    er = mask; for k=1:Opts.Erode, er = local_erode(er); end
    if nnz(er) < 50, er = mask; end                      % erosion too aggressive -> fall back

    vals = cube(er); vals = vals(isfinite(vals) & vals>0);
    vs = sort(vals); t = round(Opts.TrimPct*numel(vs));
    refTrim = mean(vs(t+1:end-t));
    switch lower(Opts.Robust)
        case 'mean',   ref = mean(vals);
        case 'median', ref = median(vals);
        otherwise,     ref = refTrim;
    end

    sMriSuvr = sMriPet;
    sMriSuvr.Cube = sMriPet.Cube ./ ref;
    sMriSuvr.Comment = [sMriPet.Comment '_suvr'];
    sMriSuvr = bst_history('add', sMriSuvr, 'suvr', sprintf( ...
        'SUVR ref=%.4g (labels %s, erode %d, %s); mask %d->%d vox', ...
        ref, mat2str(Opts.RefLabels), Opts.Erode, Opts.Robust, nnz(mask), nnz(er)));

    info = struct('RefValue',ref,'nVoxMask',nnz(mask),'nVoxEroded',nnz(er), ...
                  'RefMean',mean(vals),'RefMedian',median(vals),'RefTrim',refTrim, ...
                  'Erode',Opts.Erode,'Robust',Opts.Robust);
end

function e = local_erode(m)
    % 1-voxel 6-connectivity erosion (no Image Processing Toolbox dependency)
    e = m & circshift(m,[1 0 0]) & circshift(m,[-1 0 0]) & circshift(m,[0 1 0]) ...
          & circshift(m,[0 -1 0]) & circshift(m,[0 0 1]) & circshift(m,[0 0 -1]);
end
