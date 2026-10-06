function [sMriSuvr, info] = pet_suvr(sMriPet, sAseg, Opts)
% PET_SUVR: SUVR normalization with a robust, erosion-cleaned reference region.
%
% Upgrades the plain "mean over an ROI" rescale (mri_rescale) with the two things a reference
% region needs to be trustworthy:
%   - EROSION: peel Opts.Erode voxels off the reference mask so its boundary (which is partial-
%     volume-contaminated by neighbouring tissue/CSF) does not bias the reference downward.
%   - ROBUST estimator: a trimmed mean (default) over the eroded interior, rejecting the residual
%     outlier voxels (mis-segmentation, vessels, spill) that a plain mean is sensitive to.
% SUVR = PET / reference. Reference defaults to the cerebellar cortex (ASEG labels 8 + 47), a
% common amyloid reference (the Centiloid standard itself uses the whole cerebellum). For tau
% the inferior cerebellar GM (SUIT) is
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
%   Opts    : .Reference 'robust' (default) | 'plain'
%                 'robust': reference = 10% trimmed mean of the reference region eroded by
%                           1 voxel (positive finite voxels only). DEFAULT since this version.
%                 'plain' : reference = plain mean of all voxels of the region, no erosion:
%                           the previous behaviour (mri_rescale), for backward compatibility.
%             .RefMask (precomputed binary reference mask; overrides RefLabels/sAseg)
%             .RefLabels (default [8 47]) .Erode (voxels, 1) .Robust ('trim'|'mean'|'median')
%             .TrimPct (each-tail fraction for 'trim', 0.10). Erode/Robust/TrimPct apply
%             to 'robust' only.
%
% OUTPUTS:
%   sMriSuvr : SUVR volume (Cube ./ reference), Comment/History updated.
%   info     : struct(.RefValue,.nVoxMask,.nVoxEroded,.RefMean,.RefMedian,.RefTrim,.Erode,.Robust).
%
% SEE ALSO: pet_process, pet_pvc, mri_rescale

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

    if (nargin<3)||isempty(Opts), Opts=struct(); end
    Def=struct('Reference','robust','RefMask',[],'RefLabels',[8 47],'Erode',1,'Robust','trim','TrimPct',0.10);
    fn=fieldnames(Def); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=Def.(fn{i}); end; end

    cube = double(sMriPet.Cube(:,:,:,1));
    % Reference mask: precomputed (any region, from the caller) or built from ASEG labels.
    if ~isempty(Opts.RefMask)
        mask = logical(Opts.RefMask);
        refDesc = 'mask';
    elseif ~isempty(sAseg)
        mask = ismember(sAseg.Cube, Opts.RefLabels);
        refDesc = ['labels ' mat2str(Opts.RefLabels)];
    else
        error('pet_suvr:ref', 'Provide Opts.RefMask, or sAseg + Opts.RefLabels.');
    end
    if ~isequal(size(cube), size(mask))
        error('pet_suvr:grid', 'PET grid %s != reference mask grid %s (reslice PET to anatomy first).', ...
              mat2str(size(cube)), mat2str(size(mask)));
    end
    isPlain = strcmpi(Opts.Reference, 'plain');
    if isPlain
        % Backward compatible: plain mean over the whole region (as mri_rescale)
        er = mask;
        vals = cube(mask);
        if isempty(vals) || all(vals == 0) || ~all(isfinite(vals))
            error('pet_suvr:emptyref', 'Reference region is empty, all zero or not finite.');
        end
        Opts.Erode = 0; Opts.Robust = 'mean';
    else
        er = mask; for k=1:Opts.Erode, er = local_erode(er); end
        if nnz(er) < 50, er = mask; end                  % erosion too aggressive -> fall back
        vals = cube(er); vals = vals(isfinite(vals) & vals>0);
        if isempty(vals)
            error('pet_suvr:emptyref', 'Reference region contains no positive finite PET values.');
        end
    end
    vs = sort(vals); t = round(Opts.TrimPct*numel(vs));
    if (2*t >= numel(vs)), t = 0; end                    % too few voxels to trim
    refTrim = mean(vs(t+1:end-t));
    switch lower(Opts.Robust)
        case 'mean',   ref = mean(vals);
        case 'median', ref = median(vals);
        otherwise,     ref = refTrim;
    end

    sMriSuvr = sMriPet;
    sMriSuvr.Cube = double(sMriPet.Cube) ./ ref;         % double: an integer cube would round the SUVR
    sMriSuvr.Comment = [sMriPet.Comment '_suvr'];
    sMriSuvr = bst_history('add', sMriSuvr, 'suvr', sprintf( ...
        'SUVR ref=%.4g (%s reference, %s, erode %d, %s); mask %d->%d vox', ...
        ref, lower(Opts.Reference), refDesc, Opts.Erode, Opts.Robust, nnz(mask), nnz(er)));

    info = struct('RefValue',ref,'nVoxMask',nnz(mask),'nVoxEroded',nnz(er), ...
                  'RefMean',mean(vals),'RefMedian',median(vals),'RefTrim',refTrim, ...
                  'Erode',Opts.Erode,'Robust',Opts.Robust,'Reference',lower(Opts.Reference));
end

function e = local_erode(m)
    % 1-voxel 6-connectivity erosion (no Image Processing Toolbox dependency)
    e = m & circshift(m,[1 0 0]) & circshift(m,[-1 0 0]) & circshift(m,[0 1 0]) ...
          & circshift(m,[0 -1 0]) & circshift(m,[0 0 1]) & circshift(m,[0 0 -1]);
end
