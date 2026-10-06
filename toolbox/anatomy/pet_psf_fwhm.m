function [fwhm, src] = pet_psf_fwhm(PET, fwhmScanner)
% PET_PSF_FWHM: Effective PSF FWHM (mm) of a PET volume: scanner resolution + recorded smoothing.
%
% The scanner resolution is resolved once at import (pet_read_metadata + pet_scanner_fwhm) and
% stored as PET.PsfFwhm, so that it survives realignment, co-registration and reslicing. Any
% Gaussian smoothing applied afterwards (import "Apply smoothing", mri_realign FWHM, pet_process
% SmoothFWHM) is recorded in PET.SmoothFwhm and added in quadrature:
%
%     FWHM_eff = sqrt( FWHM_scanner^2 + sum_k FWHM_smooth(k)^2 )
%
% USAGE:  [fwhm, src] = pet_psf_fwhm(PET)
%         [fwhm, src] = pet_psf_fwhm(PET, fwhmScanner)
%
% INPUTS:
%   PET         : sMri.PET metadata struct (may be [] for volumes without metadata).
%   fwhmScanner : (optional) scanner FWHM in mm given by the user; overrides PET.PsfFwhm.
%                 The recorded smoothing is still added.
%
% OUTPUTS:
%   fwhm : effective isotropic PSF FWHM in mm.
%   src  : human-readable provenance string (for the history and the command window).
%
% SEE ALSO: pet_scanner_fwhm, pet_read_metadata, pet_gtm, pet_pvc

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

    if (nargin < 2), fwhmScanner = []; end
    % Scanner resolution: user value > value stored at import > lookup now (older files)
    if ~isempty(fwhmScanner)
        fwhm0 = double(fwhmScanner(1));
        src   = sprintf('%.1f mm (user)', fwhm0);
    elseif isstruct(PET) && isfield(PET, 'PsfFwhm') && ~isempty(PET.PsfFwhm)
        fwhm0 = double(PET.PsfFwhm(1));
        if isfield(PET, 'PsfSource') && ~isempty(PET.PsfSource)
            src = PET.PsfSource;
        else
            src = sprintf('%.1f mm (stored at import)', fwhm0);
        end
    else
        [fwhm0, src] = pet_scanner_fwhm(PET);
    end
    % Recorded smoothing, added in quadrature
    fwhm = fwhm0;
    if isstruct(PET) && isfield(PET, 'SmoothFwhm') && ~isempty(PET.SmoothFwhm)
        sm = double(PET.SmoothFwhm(:));
        sm = sm(isfinite(sm) & (sm > 0));
        if ~isempty(sm)
            fwhm = sqrt(fwhm0^2 + sum(sm.^2));
            src  = sprintf('%.2f mm = sqrt(%.1f^2 + smoothing %s^2) [scanner: %s]', ...
                           fwhm, fwhm0, mat2str(sm'), src);
        end
    end
end
