function varargout = process_source_flow( varargin )
% PROCESS_SOURCE_FLOW: Fuse a differential flow operator into an unconstrained source, producing
% a scalar flow map (divergence / curl / potential Phi / stream Psi) as a standard source file.
%
% USAGE:  OutputFiles = process_source_flow('Run', sProcess, sInputs)
%          ResultsMat = process_source_flow('Compute', ResultsMat, Quantity, Cov)
%
% The flow operators are LINEAR, so when the input is a kernel-linked source the operator is
% fused into the ImagingKernel and the output stays a kernel results file (lazy) -- which means
% the stock source analyses (PSD, Welch, average, scouts) run on the flow map directly. A
% full (ImageGridAmp) source is operated on frame-by-frame instead. Reuses the differential
% engines bst_divergence / bst_curl / bst_helmholtz_potential; geometry from the Covariant node.

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
function sProcess = GetDescription()
    sProcess.Comment     = 'Flow map (divergence/curl/potential/stream)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = '';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.options.label1.Comment = ['Fuses a differential flow operator into an unconstrained<BR>' 10 ...
        'source, producing a scalar flow map (one value per vertex).<BR>' 10 ...
        'Kernel-linked files stay kernels, so PSD / Welch / average /<BR>' 10 ...
        'scout run on the result.<BR><BR>Quantity:'];
    sProcess.options.label1.Type  = 'label';
    sProcess.options.quantity.Comment = {'<B>Divergence</B> (sources/sinks)', '<B>Curl</B> (vorticity)', ...
        '<B>Potential Phi</B> (sources/sinks, smooth)', '<B>Stream Psi</B> (vortices)'; ...
        'div', 'curl', 'phi', 'psi'};
    sProcess.options.quantity.Type  = 'radio_label';
    sProcess.options.quantity.Value = 'div';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    Quantity = sProcess.options.quantity.Value;
    tags = struct('div','flow_div', 'curl','flow_curl', 'phi','flow_phi', 'psi','flow_psi');
    fileTag = tags.(Quantity);
    nFiles = numel(sInputs);
    bst_progress('start', 'Flow map', sprintf('Computing %d files', nFiles), 0, 100);
    CovCache = containers.Map('KeyType','char', 'ValueType','any');
    for iInput = 1:nFiles
        % Keep the kernel if the source is kernel-linked (do not materialize)
        ResultsMat = in_bst_results(sInputs(iInput).FileName, 0);
        if ResultsMat.nComponents ~= 3
            bst_report('Error', sProcess, sInputs(iInput), 'Input is not an unconstrained (3-component) source model.');
            OutputFiles = {};  return;
        end
        if ~CovCache.isKey(ResultsMat.SurfaceFile)
            CovCache(ResultsMat.SurfaceFile) = tess_operators(ResultsMat.SurfaceFile, 'Covariant');
        end
        ResultsMat = Compute(ResultsMat, Quantity, CovCache(ResultsMat.SurfaceFile));
        % Save
        sStudy = bst_get('Study', sInputs(iInput).iStudy);
        OutputFiles{iInput} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['results_' fileTag]);
        bst_save(OutputFiles{iInput}, ResultsMat, 'v6');
        db_add_data(sInputs(iInput).iStudy, OutputFiles{iInput}, ResultsMat);
        bst_progress('set', round(100*iInput/nFiles));
    end
end


%% ===== COMPUTE =====
function ResultsMat = Compute(ResultsMat, Quantity, Cov)
    % Operate on the kernel (stays lazy) or the full field
    if isfield(ResultsMat,'ImagingKernel') && ~isempty(ResultsMat.ImagingKernel)
        Field = 'ImagingKernel';  X = ResultsMat.ImagingKernel;
    else
        Field = 'ImageGridAmp';   X = ResultsMat.ImageGridAmp;
    end
    switch Quantity
        case 'div',  Y = bst_divergence(X, [], 'Ambient', [], Cov);
        case 'curl', Y = bst_curl(X, [], 'Ambient', [], Cov);
        case 'phi',  [Y, ~] = bst_helmholtz_potential(X, Cov);
        case 'psi',  [~, Y] = bst_helmholtz_potential(X, Cov);
        otherwise,   error('Unknown flow quantity "%s".', Quantity);
    end
    ResultsMat.(Field)     = Y;
    ResultsMat.nComponents = 1;
    % Reset the data file initial path
    if isfield(ResultsMat,'DataFile') && ~isempty(ResultsMat.DataFile)
        ResultsMat.DataFile = file_win2unix(file_short(ResultsMat.DataFile));
    end
    labels = struct('div','flow:div', 'curl','flow:curl', 'phi','flow:Phi', 'psi','flow:Psi');
    ResultsMat.Comment = [ResultsMat.Comment ' | ' labels.(Quantity)];
    ResultsMat = bst_history('add', ResultsMat, 'flow', ['Flow map: ' Quantity]);
end
