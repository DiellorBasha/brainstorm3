function varargout = process_eigenmodes_psd( varargin )
% PROCESS_EIGENMODES_PSD: Per-eigenmode power spectrum density (Welch) of source maps.
%
% USAGE:  sProcess = process_eigenmodes_psd('GetDescription')
%       OutputFiles = process_eigenmodes_psd('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Computes the PSD of each Laplace-Beltrami mode coefficient time series by
%     substituting the mode kernel for the imaging kernel in the standard Welch
%     PSD. Output is a timefreq file: rows = modes, x = temporal frequency.
%     Requires precomputed eigenmodes on the source surface. Constrained sources.
%
% SEE ALSO: process_eigenmodes_freq, process_eigenmodes_fft, process_psd

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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Mode-frequency spectrum: PSD (Welch)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338.5;
    sProcess.Description = '';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Number of modes
    sProcess.options.nmodes.Comment = 'Number of modes: ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {300, '', 0};
    % Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Window length
    sProcess.options.win_length.Comment = 'Window length: ';
    sProcess.options.win_length.Type    = 'value';
    sProcess.options.win_length.Value   = {1, 's', []};
    % Window overlap
    sProcess.options.win_overlap.Comment = 'Window overlap ratio: ';
    sProcess.options.win_overlap.Type    = 'value';
    sProcess.options.win_overlap.Value   = {50, '%', 1};
    % Units
    sProcess.options.units.Comment = {'Physical: U<SUP>2</SUP>/Hz', 'Normalized: U<SUP>2</SUP>/Hz/s', 'Before Nov 2020', 'Units:'; ...
                                      'physical', 'normalized', 'old', ''};
    sProcess.options.units.Type    = 'radio_linelabel';
    sProcess.options.units.Value   = 'physical';
    % Save std across windows
    sProcess.options.win_std.Comment = 'Save the std across windows instead of the mean';
    sProcess.options.win_std.Type    = 'checkbox';
    sProcess.options.win_std.Value   = 0;
    % Info
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Output: timefreq file, rows = modes, x = frequency.<BR>' ...
                                           'Requires precomputed eigenmodes; constrained sources.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = process_eigenmodes_freq('Run', sProcess, sInputs, 'psd');
end
