function [fwhm, src] = pet_scanner_fwhm(PET, defaultFwhm)
% PET_SCANNER_FWHM: Estimate effective PET image resolution (PSF FWHM, mm) from metadata.
%
% BIDS PET (and the NIfTI header) carry NO explicit image-resolution / PSF field,
% so the effective in-brain FWHM is INFERRED from the scanner MODEL (a curated
% lookup of published effective resolutions) combined with any applied
% reconstruction Gaussian post-filter (ReconFilterType / ReconFilterSize). When the
% scanner cannot be identified from the metadata, a generic fallback (~6 mm, typical
% clinical PET) is returned so partial-volume correction still has a sane default.
%
% The lookup values are approximate published effective in-brain resolutions and are
% meant as sensible defaults; pass an explicit FWHM to pet_pvc to override.
%
% USAGE:  [fwhm, src] = pet_scanner_fwhm(PET)
%         [fwhm, src] = pet_scanner_fwhm(PET, defaultFwhm)
%
% INPUT:
%   PET         : sMri.PET metadata struct (see pet_read_metadata). May be [].
%   defaultFwhm : fallback FWHM in mm when the scanner is unknown (default 6).
%
% OUTPUT:
%   fwhm : estimated isotropic PSF FWHM in mm.
%   src  : human-readable provenance string (for logging / GUI).
%
% SEE ALSO: pet_read_metadata, pet_pvc
%
% Author: Diellor Basha, 2026

    if (nargin < 2) || isempty(defaultFwhm), defaultFwhm = 6; end
    fwhm = defaultFwhm;
    src  = sprintf('fallback %.1f mm (scanner not identified)', defaultFwhm);

    if isempty(PET) || ~isstruct(PET) || ~isfield(PET, 'Scanner')
        return;
    end
    model = local_str(local_getf(PET.Scanner, 'Model'));
    manuf = local_str(local_getf(PET.Scanner, 'Manufacturer'));
    hay   = lower([model ' ' manuf]);

    % Curated scanner -> intrinsic effective in-brain FWHM (mm). Order specific->general;
    % first matching key wins. Values are approximate published figures.
    tbl = { ...
        {'hrrt','high-resolution research','high resolution research'}, 2.5; ...
        {'vision quadra','biograph vision','vision'},                   3.6; ...
        {'biograph mmr','mmr'},                                         4.3; ...
        {'biograph mct','mct'},                                         4.4; ...
        {'discovery mi','dmi'},                                         4.2; ...
        {'discovery'},                                                  5.0; ...
        {'signa'},                                                      4.4; ...
        {'vereos'},                                                     4.0; ...
        {'biograph'},                                                   4.4; ...
        {'ecat'},                                                       6.0  ...
    };
    intr = []; label = '';
    for i = 1:size(tbl,1)
        keys = tbl{i,1};
        if any(cellfun(@(k) ~isempty(strfind(hay, k)), keys)) %#ok<STREMP>
            intr = tbl{i,2}; label = upper(keys{1}); break;
        end
    end
    if isempty(intr)
        return;   % keep the fallback
    end

    fwhm = intr;
    src  = sprintf('%.1f mm (%s, model lookup)', intr, label);

    % Combine with an applied reconstruction Gaussian post-filter, if present.
    ft = lower(local_str(local_getf(PET.Scanner, 'ReconFilterType')));
    fs = local_getf(PET.Scanner, 'ReconFilterSize');
    if ~isempty(fs) && isnumeric(fs) && (fs > 0) && (~isempty(strfind(ft, 'gauss'))) %#ok<STREMP>
        fwhm = sqrt(intr^2 + double(fs)^2);
        src  = sprintf('%.1f mm (%s %.1fmm + %.1fmm recon Gaussian)', fwhm, label, intr, fs);
    end
end


%% ===== helpers =====
function v = local_getf(s, f)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = []; end
end

function s = local_str(v)
    if ischar(v), s = v; elseif isempty(v), s = ''; else, s = num2str(v); end
end
