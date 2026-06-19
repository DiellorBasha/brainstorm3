function varargout = process_inverse_dirac( varargin )
% PROCESS_INVERSE_DIRAC: Source estimation in the Dirac eigenmode basis.
%
% Computes a shared imaging kernel with bst_inverse_dirac (whitened minimum-norm
% in the curvature-aware Dirac eigenbasis: amplitude / dSPM / sLORETA). The source
% space is the UNCONSTRAINED 3-vector cortical current; the Dirac transform and the
% reconstruction are handled internally, so a standard surface head model is the
% only input required (the Dirac eigenbasis is found-or-created on its surface).
%
% SEE ALSO: bst_inverse_dirac, bst_dirac, process_inverse_2018

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
    sProcess.Comment     = 'Compute sources: Dirac eigenmodes';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 327;
    sProcess.Description  = 'https://neuroimage.usc.edu/brainstorm/Tutorials/SourceEstimation';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'results', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % Options
    sProcess.options.label.Comment = ['<HTML><B>Dirac eigenmode source mapping</B><BR>' ...
        'Unconstrained whitened minimum-norm in the curvature-aware Dirac basis.<BR>' ...
        'Produces a shared imaging kernel.'];
    sProcess.options.label.Type    = 'label';
    % Inverse measure
    sProcess.options.measure.Comment = {'Current density map', 'dSPM', 'sLORETA'; 'amplitude', 'dspm2018', 'sloreta'};
    sProcess.options.measure.Type    = 'radio_label';
    sProcess.options.measure.Value   = 'dspm2018';
    % SNR
    sProcess.options.snr.Comment = 'Signal-to-noise ratio (SNR): ';
    sProcess.options.snr.Type    = 'value';
    sProcess.options.snr.Value   = {3, '', 0};
    % Modes per hemisphere
    sProcess.options.nmodes.Comment = 'Dirac modes per hemisphere: ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {400, '', 0};
    % Tau
    sProcess.options.tau.Comment = 'Dirac &tau; (intrinsic&harr;extrinsic mix): ';
    sProcess.options.tau.Type    = 'value';
    sProcess.options.tau.Value   = {0.5, '', 2};
    % Noise regularization
    sProcess.options.noisereg.Comment = 'Noise covariance regularization (fraction): ';
    sProcess.options.noisereg.Type    = 'value';
    sProcess.options.noisereg.Value   = {0.1, '', 3};
    % Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names: ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % ---- options ----
    Measure     = sProcess.options.measure.Value;
    SNR         = sProcess.options.snr.Value{1};
    nModes      = sProcess.options.nmodes.Value{1};
    Tau         = sProcess.options.tau.Value{1};
    NoiseReg    = sProcess.options.noisereg.Value{1};
    SensorTypes = sProcess.options.sensortypes.Value;

    % ---- shared kernel: one per channel study ----
    [~, iStudies] = bst_get('ChannelForStudy', unique([sInputs.iStudy]));
    iStudies = unique(iStudies);

    for iStudy = iStudies
        sStudy = bst_get('Study', iStudy);
        % --- checks ---
        if isempty(sStudy.Channel) || isempty(sStudy.Channel(1).FileName)
            bst_report('Error', sProcess, sInputs, 'No channel file in this study.'); continue;
        end
        if isempty(sStudy.HeadModel) || isempty(sStudy.iHeadModel)
            bst_report('Error', sProcess, sInputs, 'No head model selected in this study.'); continue;
        end
        if isempty(sStudy.NoiseCov) || isempty(sStudy.NoiseCov(1).FileName)
            bst_report('Error', sProcess, sInputs, 'No noise covariance in this study.'); continue;
        end
        ChannelFile = sStudy.Channel(1).FileName;
        ChannelMat  = in_bst_channel(ChannelFile);
        HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
        HeadModel = in_bst_headmodel(HeadModelFile, 0);    % full struct (incl. Dirac metadata)
        if ~strcmpi(HeadModel.HeadModelType, 'surface')
            bst_report('Error', sProcess, sInputs, 'Dirac eigenmode inverse requires a SURFACE head model.'); continue;
        end
        % Accept a pre-composed Dirac mode head model directly, or an unconstrained
        % surface head model (transformed internally). Reject scalar eigenmode models.
        isDiracHM = isfield(HeadModel,'isDiracEigenmode') && ~isempty(HeadModel.isDiracEigenmode) && HeadModel.isDiracEigenmode;
        isEigHM   = isfield(HeadModel,'isEigenmode')      && ~isempty(HeadModel.isEigenmode)      && HeadModel.isEigenmode;
        if ~isDiracHM
            if isEigHM
                bst_report('Error', sProcess, sInputs, ['Selected head model is a scalar eigenmode model. ' ...
                    'Select an "Overlapping spheres" (unconstrained) or a "Dirac eigenmode" head model.']); continue;
            end
            if mod(size(HeadModel.Gain,2), 3) ~= 0
                bst_report('Error', sProcess, sInputs, 'Head model is not an unconstrained surface model (Gain must be [nChan x 3*nVert]).'); continue;
            end
        end
        NoiseCovMat = load(file_fullpath(sStudy.NoiseCov(1).FileName));

        % --- good channels (of the requested types) ---
        ChannelFlag = ones(numel(ChannelMat.Channel), 1);
        GoodChannel = good_channel(ChannelMat.Channel, ChannelFlag, SensorTypes);
        if isempty(GoodChannel)
            bst_report('Error', sProcess, sInputs, sprintf('No good channels of type "%s".', SensorTypes)); continue;
        end

        % --- apply SSP projectors (expanded form), exactly as process_inverse_2018 ---
        if isfield(ChannelMat,'Projector') && ~isempty(ChannelMat.Projector)
            Proj = process_ssp2('BuildProjector', ChannelMat.Projector, [1 2]);
            if ~isempty(Proj)
                iGain = find(sum(isnan(HeadModel.Gain), 2) == 0);
                HeadModel.Gain(iGain,:) = Proj(iGain,iGain) * HeadModel.Gain(iGain,:);
                NoiseCovMat.NoiseCov    = Proj * NoiseCovMat.NoiseCov * Proj';
            end
        end
        % restrict to good channels
        HeadModel.Gain = HeadModel.Gain(GoodChannel, :);

        % --- options for bst_inverse_dirac ---
        OPT = struct();
        OPT.NoiseCovMat              = NoiseCovMat;
        OPT.NoiseCovMat.NoiseCov     = NoiseCovMat.NoiseCov(GoodChannel, GoodChannel);
        OPT.ChannelTypes             = {ChannelMat.Channel(GoodChannel).Type};
        OPT.NoiseMethod              = 'reg';
        OPT.NoiseReg                 = NoiseReg;
        OPT.SnrMethod                = 'fixed';
        OPT.SnrFixed                 = SNR;
        OPT.InverseMeasure           = Measure;
        OPT.nModes                   = nModes;
        OPT.Tau                      = Tau;

        % --- run the Dirac inverse (transforms + reconstructs internally) ---
        bst_progress('text', 'Computing Dirac eigenmode inverse...');
        Rd = bst_inverse_dirac(HeadModel, OPT);

        % --- FACE-based result: bst_inverse_dirac returns a per-FACE kernel [3nF x nCh];
        %     map it to the cortex VERTICES for display (each vertex = mean of incident
        %     faces, per x/y/z), mirroring how process_inverse_2018 displays face eigenmode
        %     results. The estimate stays face-native; only the display grid is the cortex. ---
        ImagingKernel = Rd.ImagingKernel;
        isFaceRes = isfield(HeadModel,'isFaceBased') && ~isempty(HeadModel.isFaceBased) && HeadModel.isFaceBased;
        if isFaceRes
            ImagingKernel = local_face_kernel_to_vertices(ImagingKernel, HeadModel.SurfaceFile);
        end

        % --- assemble the shared-kernel results structure ---
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = ImagingKernel;          % [3*nVert x nGoodChan]
        ResultsMat.nComponents   = 3;                      % unconstrained
        ResultsMat.Function      = local_function_name(Measure);
        ResultsMat.Comment       = ['Dirac: ' local_meas_comment(Measure)];
        ResultsMat.HeadModelType = 'surface';
        ResultsMat.HeadModelFile = file_short(HeadModelFile);
        ResultsMat.SurfaceFile   = file_short(HeadModel.SurfaceFile);
        ResultsMat.ChannelFlag   = ChannelFlag;
        ResultsMat.GoodChannel   = GoodChannel;
        % Persist the eigenmode-domain (amplitude) kernel + Dirac eigenvalues so the
        % Dirac mode-coefficient time-series viewer can compute c(t)=ImagingKernelMode*M
        % without re-running the inverse (see view_eigen_timeseries).
        ResultsMat.ImagingKernelMode = Rd.ImagingKernelMode;    % [nMode x nGoodChan]
        ResultsMat.Eigenvalues       = Rd.Eigenvalues;          % [nMode x 1]
        ResultsMat.ModeHemisphere    = Rd.ModeHemisphere;       % [nMode x 1]
        if isfield(Rd,'DiracEigenFile'), ResultsMat.DiracEigenFile = Rd.DiracEigenFile; end
        ResultsMat.Time          = [];
        ResultsMat.DataFile      = '';                     % shared kernel
        ResultsMat.GridLoc       = [];
        ResultsMat.GridAtlas     = [];
        ResultsMat.nAvg          = 1;
        ResultsMat.Leff          = 1;
        ResultsMat.Options       = OPT;
        ResultsMat = bst_history('add', ResultsMat, 'compute', ...
            sprintf('Dirac eigenmode source estimation (%s, SNR=%g, %d modes, tau=%g)', Measure, SNR, nModes, Tau));
        if ~isempty(sStudy.Result)
            ResultsMat.Comment = file_unique(ResultsMat.Comment, {sStudy.Result.Comment});
        end

        % --- save + register as a shared kernel ---
        OutputDir  = bst_fileparts(file_fullpath(ChannelFile));
        ResultFile = bst_process('GetNewFilename', OutputDir, 'results_DiracEig_KERNEL');
        bst_save(ResultFile, ResultsMat, 'v6');
        newResult = db_template('results');
        newResult.Comment       = ResultsMat.Comment;
        newResult.FileName      = file_short(ResultFile);
        newResult.DataFile      = '';
        newResult.isLink        = 0;
        newResult.HeadModelType = 'surface';
        sStudy.Result(end+1) = newResult;
        bst_set('Study', iStudy, sStudy);
        OutputFiles{end+1} = ResultFile; %#ok<AGROW>
    end

    % refresh the tree so the shared kernel + links appear
    if ~isempty(iStudies)
        db_links('Study', iStudies);
        panel_protocols('UpdateNode', 'Study', iStudies);
    end
end


%% ===== COMPUTE INTERACTIVE =====
% Entry point for the right-click "Compute sources: Dirac eigenmodes" menu.
% Resolves the selected tree nodes to channel studies, pops the small options
% dialog (panel_inverse_dirac), then reuses Run() to build the shared kernels.
function OutputFiles = ComputeInteractive(bstNodes) %#ok<DEFNU>
    OutputFiles = {};
    % --- resolve selected nodes to studies (mirror panel_protocols('TreeInverse')) ---
    nodeType = lower(char(bstNodes(1).getType()));
    if ismember(nodeType, {'rawdata', 'data'})
        [iStudies, ~] = tree_dependencies(bstNodes, 'data');
    else
        iStudies = tree_channel_studies(bstNodes, 'NoIntra');
    end
    if isempty(iStudies)
        return;
    elseif isequal(iStudies, -10)
        bst_error('Error in file selection.', 'Compute sources: Dirac eigenmodes', 0);
        return;
    end
    iStudies = unique(iStudies);

    % --- options dialog ---
    sOpt = gui_show_dialog('Compute sources: Dirac eigenmodes', @panel_inverse_dirac, 1, []);
    if isempty(sOpt)
        return;   % cancelled
    end

    % --- build a process structure carrying the dialog options ---
    sProcess = GetDescription();
    sProcess.options.measure.Value     = sOpt.Measure;
    sProcess.options.snr.Value{1}      = sOpt.SNR;
    sProcess.options.nmodes.Value{1}   = sOpt.nModes;
    sProcess.options.tau.Value{1}      = sOpt.Tau;
    sProcess.options.noisereg.Value{1} = sOpt.NoiseReg;
    sProcess.options.sensortypes.Value = sOpt.SensorTypes;

    % --- minimal sInputs: Run builds one shared kernel per channel study and never
    %     reads the data samples, so one entry per study (carrying iStudy) suffices ---
    sInputs = repmat(struct('iStudy', []), 1, numel(iStudies));
    for k = 1:numel(iStudies)
        sInputs(k).iStudy = iStudies(k);
    end

    % --- run the shared-kernel computation ---
    bst_progress('start', 'Compute sources', 'Computing Dirac eigenmode inverse...');
    OutputFiles = Run(sProcess, sInputs);
    bst_progress('stop');
end


%% ===== HELPERS =====
function f = local_function_name(measure)
    switch lower(measure)
        case 'amplitude', f = 'mn';
        case 'dspm2018',  f = 'dspm2018';
        case 'sloreta',   f = 'sloreta';
        otherwise,        f = measure;
    end
end

function c = local_meas_comment(measure)
    switch lower(measure)
        case 'amplitude', c = 'MN: MEG';
        case 'dspm2018',  c = 'dSPM: MEG';
        case 'sloreta',   c = 'sLORETA: MEG';
        otherwise,        c = measure;
    end
end

function Kv = local_face_kernel_to_vertices(Kf, SurfaceFile)
% Map a per-FACE unconstrained imaging kernel [3nF x nCh] to a per-VERTEX kernel
% [3nV x nCh] for display on the cortex: each vertex = mean of its incident faces,
% applied independently to the x/y/z components (rows are ordered x,y,z per element).
    Tess = in_tess_bst(SurfaceFile, 0);
    F  = double(Tess.Faces);  nV = size(Tess.Vertices,1);  nF = size(F,1);
    % vertex<-face incidence, row-normalized to the mean of incident faces
    W = sparse(F(:), repmat((1:nF)',3,1), 1, nV, nF);          % [nV x nF]
    deg = full(sum(W,2));  deg(deg==0) = 1;
    W = spdiags(1./deg, 0, nV, nV) * W;
    nCh = size(Kf,2);
    Kf3 = reshape(Kf, 3, nF, nCh);                             % [3 x nF x nCh]
    Kv  = zeros(3*nV, nCh);
    for c = 1:3
        Kv(c:3:end, :) = W * reshape(Kf3(c,:,:), nF, nCh);    % [nV x nCh]
    end
end
