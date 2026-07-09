function varargout = import_pet(varargin)
% IMPORT_PET: Import a PET volume in a subject of the Brainstorm database, and
%             read/derive its metadata.
%
% This is the single entry point for PET volumes. It reuses the full volume
% import machinery of import_mri() (loading, reorientation, PET pre-processing,
% coregistration/reslice, database registration) but presents a clean PET-only
% API, and owns the PET-specific metadata: the BIDS/NIfTI sidecar read and the
% scanner-PSF (FWHM) estimate are folded in here as sub-functions.
%
% USAGE:
%   % ---- Import a PET volume ----
%   [BstPetFile, sMri, Messages] = import_pet(iSubject, PetFile, FileFormat='ALL', isInteractive=0, isAutoAdjust=1, Comment=[])
%       - iSubject      : Index of the subject (0 = default subject)
%       - PetFile       : Full path of the PET volume (asked interactively if empty)
%       - FileFormat    : String, one of the file formats in in_mri
%       - isInteractive : If 1, importation is interactive (PET pre-processing dialog shown)
%       - isAutoAdjust  : If isInteractive=0 and isAutoAdjust=1, reslice/resample automatically
%       - Comment       : Comment of the output file
%       -> BstPetFile   : Full path to the new file if success, [] if error
%       -> sMri         : Brainstorm MRI structure (with sMri.PET metadata)
%       -> Messages     : String, messages reported by this function
%
%   % ---- Read PET metadata only (BIDS *_pet.json sidecar, NIfTI-header fallback) ----
%   PET = import_pet('ReadMetadata', PetFile, sMri)
%       -> PET : curated metadata struct stored on sMri.PET (see local_empty_pet)
%
%   % ---- Estimate effective PET PSF FWHM (mm) from scanner metadata ----
%   [fwhm, src] = import_pet('ScannerFwhm', PET, defaultFwhm=6)
%       - PET : sMri.PET metadata struct (may be [])
%       -> fwhm : estimated isotropic PSF FWHM in mm
%       -> src  : human-readable provenance string
%
% SEE ALSO: import_mri, pet_pvc, pet_gtm

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

% ===== SUB-FUNCTION DISPATCH =====
% A leading string argument selects a metadata helper; otherwise this is a volume import.
if (nargin >= 1) && ischar(varargin{1}) && ~isempty(varargin{1})
    switch (varargin{1})
        case 'ReadMetadata'
            [varargout{1:nargout}] = ReadMetadata(varargin{2:end});
            return;
        case 'ScannerFwhm'
            [varargout{1:nargout}] = ScannerFwhm(varargin{2:end});
            return;
    end
end
% Default: import a PET volume
[varargout{1:nargout}] = ImportVolume(varargin{:});
end


%% ========================================================================
%  ===== IMPORT PET VOLUME ================================================
%  ========================================================================
function [BstPetFile, sMri, Messages] = ImportVolume(iSubject, PetFile, FileFormat, isInteractive, isAutoAdjust, Comment)
    % Parse inputs
    if (nargin < 2) || isempty(PetFile)
        PetFile = [];
    end
    if (nargin < 3) || isempty(FileFormat)
        FileFormat = [];
    end
    if (nargin < 4) || isempty(isInteractive)
        isInteractive = 0;
    end
    if (nargin < 5) || isempty(isAutoAdjust)
        isAutoAdjust = 1;
    end
    if (nargin < 6)
        Comment = [];
    end
    % Tag the volume as PET so import_mri() takes its PET branch (metadata capture via
    % import_pet('ReadMetadata',...) + PET pre-processing dialog). The 'PET' prefix is
    % detected and stripped from the stored node comment by import_mri().
    if isempty(Comment)
        Comment = 'PET';
    elseif isempty(regexp(Comment, '^PET', 'once'))
        Comment = ['PET ' Comment];
    end
    % Delegate the full volume import to import_mri()
    [BstPetFile, sMri, Messages] = import_mri(iSubject, PetFile, FileFormat, isInteractive, isAutoAdjust, Comment);
end


%% ========================================================================
%  ===== READ PET METADATA ================================================
%  ========================================================================
% Build the curated PET metadata struct (sMri.PET) for a PET volume.
%
% Unlike a structural MRI, a PET volume is uninterpretable without its tracer,
% frame/injection timing and scanner/recon context. The primary source is the
% BIDS JSON sidecar next to PetFile (*_pet.json); if absent, a partial struct is
% derived from the NIfTI header / loaded cube (sMri). Pure read/parse: no DB side effects.
function PET = ReadMetadata(PetFile, sMri)
    if (nargin < 2)
        sMri = [];
    end
    PET = local_empty_pet();
    % ---- primary source: BIDS JSON sidecar ----
    jsonFile = local_find_sidecar(PetFile);
    if ~isempty(jsonFile)
        try
            J = bst_jsondecode(jsonFile);
        catch
            J = [];
        end
        if ~isempty(J)
            PET = local_from_json(PET, J);
            PET.Source = 'bids-json';
        end
    end
    % ---- fallback: NIfTI header / loaded cube ----
    if strcmp(PET.Source, 'none')
        PET = local_from_header(PET, sMri);
    end
    % ---- derived frame fields ----
    PET = local_derive_frames(PET);
end

% ===== empty schema =====
function PET = local_empty_pet()
    PET = struct();
    PET.Source    = 'none';
    PET.Tracer    = struct('Name','', 'Radionuclide','', 'Units','');
    PET.Injection = struct('InjectedRadioactivity',[], 'InjectedRadioactivityUnits','', ...
                           'Mode','', 'TimeZero','', 'InjectionStart',[], 'ScanStart',[]);
    PET.Frames    = struct('TimesStart',[], 'Duration',[], 'N',0, 'MidTimes',[], 'CoverageMinPI',[]);
    PET.Decay     = struct('ImageDecayCorrected',[], 'ImageDecayCorrectionTime',[]);
    PET.Scanner   = struct('Manufacturer','', 'Model','', 'ReconMethod','', ...
                           'ReconFilterType','', 'ReconFilterSize',[], 'AttenuationCorrection','');
    PET.Json      = [];
end

% ===== locate the BIDS JSON sidecar =====
function jsonFile = local_find_sidecar(PetFile)
    jsonFile = '';
    if isempty(PetFile) || iscell(PetFile)
        return;
    end
    [fPath, fBase, fExt] = bst_fileparts(PetFile);
    if strcmpi(fExt, '.gz')          % *.nii.gz -> strip the .nii too
        [~, fBase] = bst_fileparts(fBase);
    end
    if strncmp(fBase, '._', 2)       % macOS AppleDouble companion
        return;
    end
    cand = bst_fullfile(fPath, [fBase '.json']);
    if file_exist(cand)
        jsonFile = cand;
    end
end

% ===== map BIDS JSON -> curated schema =====
function PET = local_from_json(PET, J)
    PET.Json = J;
    PET.Tracer.Name         = local_get(J, 'TracerName', '');
    PET.Tracer.Radionuclide = local_get(J, 'TracerRadionuclide', '');
    PET.Tracer.Units        = local_get(J, 'Units', '');

    PET.Injection.InjectedRadioactivity      = local_get(J, 'InjectedRadioactivity', []);
    PET.Injection.InjectedRadioactivityUnits = local_get(J, 'InjectedRadioactivityUnits', '');
    PET.Injection.Mode           = local_get(J, 'ModeOfAdministration', '');
    PET.Injection.TimeZero       = local_get(J, 'TimeZero', '');
    PET.Injection.InjectionStart = local_get(J, 'InjectionStart', []);
    PET.Injection.ScanStart      = local_get(J, 'ScanStart', []);

    PET.Frames.TimesStart = local_row(local_get(J, 'FrameTimesStart', []));
    PET.Frames.Duration   = local_row(local_get(J, 'FrameDuration', []));

    PET.Decay.ImageDecayCorrected      = local_get(J, 'ImageDecayCorrected', []);
    PET.Decay.ImageDecayCorrectionTime = local_get(J, 'ImageDecayCorrectionTime', []);

    PET.Scanner.Manufacturer          = local_get(J, 'Manufacturer', '');
    PET.Scanner.Model                 = local_get(J, 'ManufacturersModelName', '');
    PET.Scanner.ReconMethod           = local_get(J, 'ReconMethodName', '');
    PET.Scanner.ReconFilterType       = local_get(J, 'ReconFilterType', '');
    PET.Scanner.ReconFilterSize       = local_get(J, 'ReconFilterSize', []);
    PET.Scanner.AttenuationCorrection = local_get(J, 'AttenuationCorrection', '');
end

% ===== NIfTI-header / cube fallback =====
function PET = local_from_header(PET, sMri)
    if isempty(sMri)
        return;   % Source stays 'none'
    end
    % Frame count from the loaded cube (most reliable indicator of N frames)
    if isfield(sMri, 'Cube') && ~isempty(sMri.Cube)
        PET.Frames.N = size(sMri.Cube, 4);
    end
    % Units, if the header carries them
    if isfield(sMri, 'Header') && ~isempty(sMri.Header) && isfield(sMri.Header, 'nifti') ...
            && isfield(sMri.Header.nifti, 'intent_name') && ~isempty(sMri.Header.nifti.intent_name)
        PET.Tracer.Units = sMri.Header.nifti.intent_name;
    end
    PET.Source = 'nifti-header';
end

% ===== derived frame fields =====
function PET = local_derive_frames(PET)
    ts = PET.Frames.TimesStart;
    du = PET.Frames.Duration;
    if isempty(ts)
        return;
    end
    PET.Frames.N = numel(ts);
    if ~isempty(du) && (numel(du) == numel(ts))
        PET.Frames.MidTimes = ts + du / 2;
    end
    inj = PET.Injection.InjectionStart;
    if ~isempty(inj) && ~isempty(du) && (numel(du) == numel(ts))
        t0 = (min(ts)        - inj) / 60;    % minutes post-injection
        t1 = (max(ts + du)   - inj) / 60;
        PET.Frames.CoverageMinPI = [t0, t1];
    end
end

% ===== getters =====
function v = local_get(J, name, def)
% Field value with default; BIDS 'n/a' sentinels map to the default.
    if isfield(J, name) && ~isempty(J.(name))
        v = J.(name);
        if ischar(v) && strcmpi(strtrim(v), 'n/a')
            v = def;
        end
    else
        v = def;
    end
end

function v = local_row(x)
% Coerce a JSON numeric array (column vector from jsondecode, or cell from the
% Brainstorm fallback decoder) into a double row vector.
    if iscell(x)
        x = cell2mat(x);
    end
    if isempty(x)
        v = [];
    else
        v = double(x(:)');
    end
end


%% ========================================================================
%  ===== SCANNER PSF FWHM =================================================
%  ========================================================================
% Estimate effective PET image resolution (PSF FWHM, mm) from metadata.
%
% BIDS PET (and the NIfTI header) carry NO explicit image-resolution / PSF field,
% so the effective in-brain FWHM is INFERRED from the scanner MODEL (a curated
% lookup of published effective resolutions) combined with any applied
% reconstruction Gaussian post-filter. When the scanner cannot be identified, a
% generic fallback (~6 mm, typical clinical PET) is returned. Pass an explicit
% FWHM to pet_pvc to override.
function [fwhm, src] = ScannerFwhm(PET, defaultFwhm)
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

% ===== scanner-struct helpers =====
function v = local_getf(s, f)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = []; end
end

function s = local_str(v)
    if ischar(v), s = v; elseif isempty(v), s = ''; else, s = num2str(v); end
end
