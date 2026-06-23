function varargout = process_eigen( varargin )
% PROCESS_EIGEN: Pure orchestrator for eigen-domain analysis (the spatial-spectral twin of
% PROCESS_TIMEFREQ). It is NOT a registered process (no GetDescription => never shown in the
% pipeline panel). Sibling processes (process_eigenspectrum, and future process_eigenfilter/
% process_eigenwavelet) delegate their Run to this function, exactly as process_psd/
% process_hilbert delegate to process_timefreq('Run', ...). The eigen METHOD is inferred from
% the calling process name; panel options are translated to bst_eigen OPTIONS; bst_eigen does
% the read -> resolve eigen_ basis -> compute -> save.

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


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % Method from the calling process name (mirrors process_timefreq.m:91-100)
    switch func2str(sProcess.Function)
        case 'process_eigenspectrum', Method = 'spectrum';
        case 'process_eigenfilter',   Method = 'filter';
        case 'process_eigenwavelet',  Method = 'wavelet';
        otherwise
            bst_report('Error', sProcess, sInputs, 'Unsupported eigen process.');
            return;
    end
    % Build engine options: defaults, then overlay the panel options.
    OPTIONS = bst_eigen();
    OPTIONS.Method       = Method;
    OPTIONS.iTargetStudy = [];   % each output node lands in its own input's study
    % Variant (radio_linelabel: Value is the key string), default Laplace-Beltrami
    if isfield(sProcess.options, 'variant') && ~isempty(sProcess.options.variant.Value)
        OPTIONS.Variant = sProcess.options.variant.Value;
    end
    % Measure (power/magnitude)
    if isfield(sProcess.options, 'measure') && ~isempty(sProcess.options.measure.Value)
        OPTIONS.Measure = sProcess.options.measure.Value;
    end
    % Time window
    if isfield(sProcess.options, 'timewindow') && ~isempty(sProcess.options.timewindow.Value) ...
            && iscell(sProcess.options.timewindow.Value)
        OPTIONS.TimeWindow = sProcess.options.timewindow.Value{1};
    end
    % Welch-style windowing (window length empty => single full window)
    if isfield(sProcess.options, 'win_length') && ~isempty(sProcess.options.win_length.Value) ...
            && iscell(sProcess.options.win_length.Value)
        OPTIONS.WinLength  = sProcess.options.win_length.Value{1};
        OPTIONS.WinOverlap = sProcess.options.win_overlap.Value{1};
    end
    % Aggregate across windows
    if isfield(sProcess.options, 'win_std') && ~isempty(sProcess.options.win_std.Value)
        if sProcess.options.win_std.Value
            OPTIONS.WinFunc = 'mean+std';
        else
            OPTIONS.WinFunc = 'mean';
        end
    end
    % Run the engine per input so one failure does not abort the batch.
    for iIn = 1:numel(sInputs)
        try
            [out, Messages, isError] = bst_eigen(sInputs(iIn).FileName, OPTIONS);
            if isError
                bst_report('Error', sProcess, sInputs(iIn), Messages);
            else
                if ~isempty(Messages)
                    bst_report('Warning', sProcess, sInputs(iIn), Messages);
                end
                OutputFiles = [OutputFiles, out]; %#ok<AGROW>
            end
        catch ME
            bst_report('Error', sProcess, sInputs(iIn), ME.message);
        end
    end
end
