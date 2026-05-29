function varargout = process_eigenmodes_freq( varargin )
% PROCESS_EIGENMODES_FREQ: Shared engine for the mode-frequency spectrum (PSD/FFT).
%
% Not listed in the process menu (no GetDescription). Driven by the sibling
% processes process_eigenmodes_psd / process_eigenmodes_fft.
%
% USAGE:  OutputFiles = process_eigenmodes_freq('Run', sProcess, sInputs, Method)
%         [TF, Freqs, RowNames, Msg] = process_eigenmodes_freq('Compute', ...
%               F, sfreq, ModeKernel, BadSegments, Method, WinLength, WinOverlap, WinFunc, PowerUnits, Values)
%
% DESCRIPTION:
%     Computes the per-eigenmode temporal-frequency spectrum of a source map by
%     substituting the mode kernel (Phi'*M*ImagingKernel) for the imaging kernel
%     in the standard bst_psd call. Because bst_psd applies the kernel after the
%     FFT (valid by linearity) and skips bad-segment windows before the kernel,
%     the result is identical to a source PSD/FFT with Phi'*M applied.
%
% SEE ALSO: bst_eigenmodes_modekernel, bst_psd, process_eigenmodes_psd, process_eigenmodes_fft

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


%% ===== COMPUTE (pure spectral core) =====
function [TF, Freqs, RowNames, Msg] = Compute(F, sfreq, ModeKernel, BadSegments, Method, WinLength, WinOverlap, WinFunc, PowerUnits, Values) %#ok<DEFNU>
    nModes = size(ModeKernel, 1);
    switch lower(Method)
        case 'psd'
            [TF, Freqs, ~, Msg] = bst_psd(F, sfreq, WinLength, WinOverlap, BadSegments, ModeKernel, WinFunc, PowerUnits);
        case 'fft'
            [TF, Freqs, ~, Msg] = bst_psd(F, sfreq, [], 0, BadSegments, ModeKernel, [], PowerUnits);
        otherwise
            error('Unsupported method: %s', Method);
    end
    % Row labels: mode index + eigenvalue
    RowNames = cell(nModes, 1);
    for k = 1:nModes
        RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, Values(k));
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs, Method) %#ok<DEFNU>
    OutputFiles = {};

    % ===== OPTIONS =====
    nModesReq = sProcess.options.nmodes.Value{1};
    PowerUnits = 'physical';
    if isfield(sProcess.options, 'units') && ~isempty(sProcess.options.units.Value)
        PowerUnits = sProcess.options.units.Value;
    end
    TimeWindow = [];
    if isfield(sProcess.options, 'timewindow') && ~isempty(sProcess.options.timewindow.Value) && iscell(sProcess.options.timewindow.Value)
        TimeWindow = sProcess.options.timewindow.Value{1};
    end
    if strcmpi(Method, 'psd')
        WinLength  = sProcess.options.win_length.Value{1};
        WinOverlap = sProcess.options.win_overlap.Value{1};
        switch (sProcess.options.win_std.Value)
            case {0, 'mean'}, WinFunc = 'mean';
            case {1, 'std'},  WinFunc = 'std';
            otherwise,        WinFunc = 'mean';
        end
    else
        WinLength = []; WinOverlap = 0; WinFunc = 'mean';
    end

    bst_progress('start', 'Mode-frequency spectrum', 'Computing...');

    % ===== PER INPUT FILE =====
    for iInput = 1:length(sInputs)
        sInput = sInputs(iInput);

        % --- Load metadata only (no kernel expansion) ---
        ResultsMat = in_bst_results(sInput.FileName, 0, ...
            'ImagingKernel', 'ImageGridAmp', 'GoodChannel', 'nComponents', ...
            'DataFile', 'Time', 'SurfaceFile', 'HeadModelType', 'Atlas');

        % --- Validate ---
        if isfield(ResultsMat, 'HeadModelType') && ~isempty(ResultsMat.HeadModelType) && ~strcmpi(ResultsMat.HeadModelType, 'surface')
            bst_report('Error', sProcess, sInput, 'Only surface source models are supported.'); continue;
        end
        if isfield(ResultsMat, 'Atlas') && ~isempty(ResultsMat.Atlas)
            bst_report('Error', sProcess, sInput, 'Atlas-based source models are not supported.'); continue;
        end
        if ResultsMat.nComponents ~= 1
            bst_report('Error', sProcess, sInput, 'Constrained source orientation required (nComponents == 1).'); continue;
        end
        SurfaceFile = ResultsMat.SurfaceFile;
        if isempty(SurfaceFile)
            bst_report('Error', sProcess, sInput, 'No surface file associated with the source map.'); continue;
        end

        % --- Eigenmodes ---
        [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
        if ~isComputed
            bst_report('Error', sProcess, sInput, ['No eigenmodes on: ' SurfaceFile '. Run "Compute eigenmodes" first.']); continue;
        end
        nV_eigen = size(Eig.Vectors, 1);

        % --- Mass matrix (reuse if saved) ---
        if isfield(Eig, 'MassMatrix') && ~isempty(Eig.MassMatrix)
            M = Eig.MassMatrix;
        else
            sSurf = in_tess_bst(SurfaceFile, 0);
            [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
        end

        % --- Mode count ---
        nModes = max(1, min(nModesReq, size(Eig.Vectors, 2)));
        if nModes < nModesReq
            bst_report('Warning', sProcess, sInput, sprintf('Requested %d modes but only %d available; using %d.', nModesReq, nModes, nModes));
        end

        % --- Read sensor/source data + bad segments (mirrors bst_timefreq) ---
        [F, TimeVector, BadSegments, errMsg] = ReadInput(sInput, ResultsMat, TimeWindow, Method);
        if ~isempty(errMsg)
            bst_report('Error', sProcess, sInput, errMsg); continue;
        end
        sfreq = 1 ./ (TimeVector(2) - TimeVector(1));

        % --- Build the kernel passed to bst_psd ---
        if isempty(ResultsMat.ImagingKernel)
            % Full result: F are sources; check vertex match; kernel = projector
            if size(F,1) ~= nV_eigen
                bst_report('Error', sProcess, sInput, sprintf('Vertex mismatch: sources %d, eigenmodes %d.', size(F,1), nV_eigen)); continue;
            end
            ModeKernel = bst_eigenmodes_modekernel(Eig, M, [], nModes);            % [nModes x nV]
        else
            % Kernel result: F are sensors; kernel = Phi'*M*ImagingKernel
            ModeKernel = bst_eigenmodes_modekernel(Eig, M, ResultsMat.ImagingKernel, nModes);  % [nModes x nChan]
        end

        % --- Compute the spectrum (kernel-swap) ---
        [TF, Freqs, RowNames, Msg] = Compute(F, sfreq, ModeKernel, BadSegments, Method, WinLength, WinOverlap, WinFunc, PowerUnits, Eig.Values);
        if isempty(TF)
            bst_report('Error', sProcess, sInput, ['Spectral computation failed: ' Msg]); continue;
        end

        % --- Build timefreq file ---
        FileMat = db_template('timefreqmat');
        FileMat.TF          = TF;
        FileMat.Time        = TimeVector([1 end]);
        FileMat.Freqs       = Freqs;
        FileMat.RowNames    = RowNames;
        FileMat.Measure     = 'power';
        FileMat.Method      = Method;
        FileMat.DataType    = 'matrix';
        FileMat.SurfaceFile = SurfaceFile;
        FileMat.DataFile    = '';
        FileMat.nAvg        = 1;
        FileMat.Leff        = 1;
        switch lower(Method)
            case 'psd', FileMat.Comment = sprintf('Mode PSD [%d modes]', nModes);
            case 'fft', FileMat.Comment = sprintf('Mode FFT [%d modes]', nModes);
        end
        FileMat.Options.Method     = Method;
        FileMat.Options.Measure    = 'power';
        FileMat.Options.PowerUnits = PowerUnits;
        if strcmpi(Method, 'psd')
            FileMat.Options.WindowFunction = WinFunc;
        end
        FileMat = bst_history('add', FileMat, 'compute', sprintf('Mode-frequency spectrum (%s) of: %s', Method, sInput.FileName));
        FileMat = bst_history('add', FileMat, 'compute', sprintf('%d modes, %s mass, eigenvalue range [%.2f, %.2f]', nModes, Eig.MassType, Eig.Values(1), Eig.Values(nModes)));

        % --- Save ---
        [sStudy, iStudy] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
        OutputFile = bst_process('GetNewFilename', StudyDir, ['timefreq_' Method]);
        bst_save(OutputFile, FileMat, 'v6');
        db_add_data(iStudy, OutputFile, FileMat);
        OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>
    end

    bst_progress('stop');
end


%% ===== READ INPUT (sensor/source data + bad segments) =====
% Mirrors the results branch of bst_timefreq (bst_timefreq.m:346-440) so that
% bad-segment handling and the time-window read are identical to a source PSD.
function [F, TimeVector, BadSegments, errMsg] = ReadInput(sInput, ResultsMat, TimeWindow, Method)
    F = []; TimeVector = []; BadSegments = []; errMsg = '';

    if ~isempty(ResultsMat.ImagingKernel) && isempty(ResultsMat.ImageGridAmp)
        % ----- Kernel result: process the recordings file -----
        DataFile = ResultsMat.DataFile;
        if isempty(DataFile)
            errMsg = 'Kernel result requires an associated data file.'; return;
        end
        sMat = in_bst_data(DataFile);
        if isstruct(sMat.F)
            % Raw recordings: bounded windowed read
            sFile = sMat.F;
            ChannelFile = bst_get('ChannelFileForStudy', sInput.FileName);
            if isempty(ChannelFile)
                errMsg = 'No channel definition available for this file.'; return;
            end
            ChannelMat = in_bst_channel(ChannelFile);
            if (length(sFile.epochs) > 1)
                errMsg = 'Files with epochs are not supported by this process.'; return;
            end
            [F, TimeVector, BadSegments] = ReadRawRecordings(sFile, sMat.Time, ChannelMat, TimeWindow, Method);
        else
            % Imported recordings
            F = sMat.F;
            TimeVector = sMat.Time;
            sMat.events = sMat.Events;
            sMat.prop.sfreq = 1 ./ (sMat.Time(2) - sMat.Time(1));
            isChannelEvtBad = 0;
            BadSegments = panel_record('GetBadSegments', sMat, isChannelEvtBad) - sMat.prop.sfreq * sMat.Time(1) + 1;
        end
        % Restrict to this result's good channels
        F = F(ResultsMat.GoodChannel, :);
    else
        % ----- Full result: sources themselves -----
        F = ResultsMat.ImageGridAmp;
        TimeVector = ResultsMat.Time;
        BadSegments = [];
    end

    % ----- Keep only the requested time window (mirrors bst_timefreq.m:425-440) -----
    if ~isempty(TimeWindow)
        iTime = bst_closest(TimeWindow, TimeVector);
        if (iTime(1) == iTime(2))
            errMsg = 'Selected time window is not valid for the input file.'; return;
        end
        iTime = iTime(1):iTime(2);
        TimeVector = TimeVector(iTime);
        F = F(:, iTime);
    end
end


%% ===== READ RAW RECORDINGS (bounded window) =====
% Copy of bst_timefreq's private ReadRawRecordings (bst_timefreq.m:960), with
% TimeWindow passed in explicitly. PSD bad-segment handling preserved.
function [F, TimeVector, BadSegments] = ReadRawRecordings(sFile, TimeVector, ChannelMat, TimeWindow, Method)
    ImportOptions = db_template('ImportOptions');
    ImportOptions.ImportMode     = 'Time';
    ImportOptions.Resample       = 0;
    ImportOptions.UseCtfComp     = 1;
    ImportOptions.UseSsp         = 1;
    ImportOptions.RemoveBaseline = 'no';
    ImportOptions.DisplayMessages = 0;
    % Samples to read
    if ~isempty(TimeWindow)
        SamplesBounds = round(sFile.prop.times(1) .* sFile.prop.sfreq) + bst_closest(TimeWindow, TimeVector) - 1;
    else
        SamplesBounds = round(sFile.prop.times .* sFile.prop.sfreq);
    end
    % Read data
    [F, TimeVector] = in_fread(sFile, ChannelMat, 1, SamplesBounds, [], ImportOptions);
    % Bad segments: only used by PSD (Welch averaging). FFT is a single
    % full-recording window, so bad segments are not applied (mirrors
    % bst_timefreq, which gates GetBadSegments on the 'psd' method).
    if strcmpi(Method, 'psd')
        isChannelEvtBad = 0;
        BadSegments = panel_record('GetBadSegments', sFile, isChannelEvtBad);
        if ~isempty(BadSegments)
            BadSegments = BadSegments - SamplesBounds(1) + 1;
        end
    else
        BadSegments = [];
    end
end
