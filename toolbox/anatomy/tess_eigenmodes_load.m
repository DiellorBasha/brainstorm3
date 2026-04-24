function [Eigenmodes, isComputed] = tess_eigenmodes_load(SurfaceFile)
% TESS_EIGENMODES_LOAD: Load precomputed eigenmodes from a Brainstorm surface file.
%
% USAGE:  [Eigenmodes, isComputed] = tess_eigenmodes_load(SurfaceFile)
%
% DESCRIPTION:
%     Loads the Eigenmodes field from a Brainstorm surface file. Returns
%     empty and false if eigenmodes have not been computed for this surface.
%
% INPUTS:
%     SurfaceFile : Relative path to the surface file in the Brainstorm
%                   database, or full absolute path.
%
% OUTPUTS:
%     Eigenmodes  : Structure with fields .Vectors, .Values, .nModes, etc.
%                   Returns [] if not computed.
%     isComputed  : Logical, true if eigenmodes exist in the file.
%
% SEE ALSO: tess_eigenmodes, tess_eigenmodes_save

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

Eigenmodes = [];
isComputed = false;

%% ===== RESOLVE PATH =====
if ~isempty(strfind(SurfaceFile, filesep)) && exist(SurfaceFile, 'file') %#ok<STREMP>
    SurfaceFileFull = SurfaceFile;
else
    SurfaceFileFull = file_fullpath(SurfaceFile);
end

if ~exist(SurfaceFileFull, 'file')
    error('Surface file not found: %s', SurfaceFileFull);
end

%% ===== CHECK FOR EIGENMODES =====
% Use matfile for fast check without loading the entire surface
try
    varInfo = who('-file', SurfaceFileFull);
    if ~ismember('Eigenmodes', varInfo)
        return;
    end

    % Load only the Eigenmodes field
    S = load(SurfaceFileFull, 'Eigenmodes');
    if ~isfield(S, 'Eigenmodes') || isempty(S.Eigenmodes)
        return;
    end

    Eigenmodes = S.Eigenmodes;

    % Convert Vectors back to double for computation
    if isa(Eigenmodes.Vectors, 'single')
        Eigenmodes.Vectors = double(Eigenmodes.Vectors);
    end

    isComputed = true;
catch
    % File doesn't have Eigenmodes or is corrupted
    return;
end

end
