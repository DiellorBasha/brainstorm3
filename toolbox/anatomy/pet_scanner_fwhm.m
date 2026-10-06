function [fwhm, src] = pet_scanner_fwhm(PET, defaultFwhm)
% PET_SCANNER_FWHM: Estimate effective PET image resolution (PSF FWHM, mm) from metadata.
%
% BIDS PET (and the NIfTI header) carry NO explicit image-resolution / PSF field,
% so the FWHM is INFERRED from the scanner MODEL (a curated lookup of the published
% NEMA spatial resolution, matched exactly on the model name) combined with any applied
% reconstruction Gaussian post-filter (ReconFilterType / ReconFilterSize). When the
% scanner cannot be identified from the metadata, a generic fallback (~6 mm, typical
% clinical PET) is returned so partial-volume correction still has a sane default.
%
% The lookup values are published NEMA resolutions (references in the table below) and
% are meant as sensible defaults; pass an explicit FWHM to pet_pvc to override.
% Smoothing applied after reconstruction is NOT included here: see pet_psf_fwhm.
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
% SEE ALSO: pet_psf_fwhm, pet_read_metadata, pet_pvc

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

    if (nargin < 2) || isempty(defaultFwhm), defaultFwhm = 6; end
    fwhm = defaultFwhm;
    src  = sprintf('fallback %.1f mm (scanner not identified)', defaultFwhm);

    if isempty(PET) || ~isstruct(PET) || ~isfield(PET, 'Scanner')
        return;
    end
    model = local_str(local_getf(PET.Scanner, 'Model'));
    % Exact matching on the normalized model name (lower case, letters and digits only), so that
    % a short key cannot match an unrelated scanner (e.g. "mct" inside another model name).
    key = lower(regexprep(model, '[^a-zA-Z0-9]', ''));

    % Curated scanner -> intrinsic resolution (mm): NEMA NU-2 transaxial spatial resolution
    % (FWHM at 1 cm from the centre of the field of view), from the published performance
    % evaluation of each system. The effective in-brain resolution of a reconstructed image
    % can be worse (reconstruction, motion); pass an explicit FWHM to override.
    %  - HRRT: de Jong et al. (2007), Phys Med Biol 52:1505-1526
    %  - Biograph Vision: van Sluis et al. (2019), J Nucl Med 60:1031-1036
    %  - Biograph mMR: Delso et al. (2011), J Nucl Med 52:1914-1922
    %  - Biograph mCT: Jakoby et al. (2011), Phys Med Biol 56:2375-2389
    %  - Discovery MI: Hsu et al. (2017), J Nucl Med 58:1511-1518
    %  - SIGNA PET/MR: Grant et al. (2016), Med Phys 43:2334-2343
    %  - Vereos: Zhang et al. (2018), EJNMMI Res 8:97
    tbl = { ...
        'HRRT',            {'hrrt', 'ecathrrt', 'highresolutionresearchtomograph'},                     2.5; ...
        'Biograph Vision', {'biographvision', 'biographvision600', 'biographvisionquadra'},             3.6; ...
        'Biograph mMR',    {'biographmmr', 'mmr'},                                                      4.3; ...
        'Biograph mCT',    {'biographmct', 'biograph128mct', 'biograph64mct', 'biographmctflow', 'mct'}, 4.4; ...
        'Discovery MI',    {'discoverymi', 'discoverymidr'},                                            4.2; ...
        'SIGNA PET/MR',    {'signapetmr', 'signa'},                                                     4.4; ...
        'Vereos',          {'vereos', 'vereospetct'},                                                   4.0  ...
    };
    intr = []; label = '';
    if ~isempty(key)
        for i = 1:size(tbl,1)
            if ismember(key, tbl{i,2})
                intr = tbl{i,3}; label = tbl{i,1}; break;
            end
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
    if isnumeric(fs) && isscalar(fs) && (fs > 0) && (~isempty(strfind(ft, 'gauss'))) %#ok<STREMP>
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
