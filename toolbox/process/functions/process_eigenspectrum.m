function varargout = process_eigenspectrum( varargin )
% PROCESS_EIGENSPECTRUM: Windowed eigenspectrum of surface-mapped sources (the spatial-spectral
% analogue of PROCESS_PSD / Welch). Projects each source map onto an operator eigenbasis and
% returns mode power vs spatial frequency (Freqs = sqrt(lambda)) as a timefreq_ node. The eigen_
% basis is resolved implicitly from the input's SurfaceFile; the operator family is chosen by the
% Variant option. Thin wrapper: delegates Run to the pure orchestrator process_eigen, exactly as
% process_psd delegates to process_timefreq.

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
    sProcess.Comment     = 'Eigenspectrum (spatial PSD)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Frequency';
    sProcess.Index       = 484;
    sProcess.Description = '';
    % Source maps only (eigenspectrum needs a surface basis)
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Option: operator family (eigen_ basis resolved implicitly from the surface)
    sProcess.options.variant.Comment = {'Laplace-Beltrami', 'Connection Laplacian', 'Dirac', 'Dirac-Face', 'Hodge-Face', 'Operator:'; ...
        'Laplace-Beltrami', 'Connection Laplacian', 'Dirac', 'Dirac-Face', 'Hodge-Face', ''};
    sProcess.options.variant.Type    = 'radio_linelabel';
    sProcess.options.variant.Value   = 'Laplace-Beltrami';
    % Option: spectrum measure
    sProcess.options.measure.Comment = {'Power', 'Magnitude', 'Measure:'; 'power', 'magnitude', ''};
    sProcess.options.measure.Type    = 'radio_linelabel';
    sProcess.options.measure.Value   = 'power';
    % Option: time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Option: Welch window length (empty => single full window)
    sProcess.options.win_length.Comment = 'Window length: ';
    sProcess.options.win_length.Type    = 'value';
    sProcess.options.win_length.Value   = {1, 's', []};
    % Option: Welch window overlap
    sProcess.options.win_overlap.Comment = 'Window overlap ratio: ';
    sProcess.options.win_overlap.Type    = 'value';
    sProcess.options.win_overlap.Value   = {50, '%', 1};
    % Option: aggregate (mean vs mean+std across windows)
    sProcess.options.win_std.Comment = '<HTML><FONT color="#a0a0a0">Also save the std across windows</FONT>';
    sProcess.options.win_std.Type    = 'checkbox';
    sProcess.options.win_std.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Delegate to the shared eigen orchestrator (method inferred from this process name).
    OutputFiles = process_eigen('Run', sProcess, sInputs);
end
