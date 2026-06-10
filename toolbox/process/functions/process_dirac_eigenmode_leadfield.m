function varargout = process_dirac_eigenmode_leadfield( varargin )
% PROCESS_DIRAC_EIGENMODE_LEADFIELD: Compose the unconstrained leadfield into the Dirac eigenbasis.
%
% Produces a headmodel_dirac_eigenmode_*.mat node whose Gain [nCh x 2K] expresses
% the UNCONSTRAINED leadfield in the Dirac (curvature-aware vector) eigenbasis
% (Gain = Psi' * B * Phi_D). Requires a surface unconstrained head model; computes
% the Dirac eigenbasis on the cortex (tess_dirac_eigenmodes) if absent.
%
% Authors: Diellor Basha, 2026
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Compute Dirac eigenmode leadfield';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 340;
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    sProcess.options.nmodes.Comment = 'Number of Dirac eigenmodes per hemisphere (0 = default 400): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};
    sProcess.options.tau.Comment    = 'Curvature weighting tau [0-1]: ';
    sProcess.options.tau.Type       = 'value';
    sProcess.options.tau.Value      = {0.5, '', 3};
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Expresses the UNCONSTRAINED leadfield in the ' ...
        'Dirac (curvature-aware vector) eigenbasis.<BR>Requires a surface head model; computes the Dirac ' ...
        'eigenbasis on the cortex if absent.</FONT>'];
    sProcess.options.label_info.Type = 'label';
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    n = sProcess.options.nmodes.Value{1};
    if n > 0; Comment = sprintf('Dirac eigenmode leadfield (%d modes)', n);
    else;     Comment = 'Dirac eigenmode leadfield (default modes)'; end
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    nModes = sProcess.options.nmodes.Value{1};
    Tau    = sProcess.options.tau.Value{1};
    Kbasis = nModes; if Kbasis <= 0, Kbasis = 400; end

    [sStudy, iStudy] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.HeadModel)
        bst_report('Error', sProcess, sInputs, 'No head model available for this study.'); return;
    end

    % Choose a BASE (non-eigenmode) UNCONSTRAINED surface head model.
    iBase = []; HeadModel = [];
    candidates = [sStudy.iHeadModel, setdiff(1:numel(sStudy.HeadModel), sStudy.iHeadModel)];
    for ic = candidates
        if ic < 1 || ic > numel(sStudy.HeadModel); continue; end
        try, hmC = in_bst_headmodel(sStudy.HeadModel(ic).FileName, 0); catch, continue; end
        isEig  = isfield(hmC,'isEigenmode')      && ~isempty(hmC.isEigenmode)      && hmC.isEigenmode;
        isDEig = isfield(hmC,'isDiracEigenmode') && ~isempty(hmC.isDiracEigenmode) && hmC.isDiracEigenmode;
        if ~isEig && ~isDEig && strcmpi(hmC.HeadModelType, 'surface')
            iBase = ic; HeadModel = hmC; break;
        end
    end
    if isempty(iBase)
        bst_report('Error', sProcess, sInputs, 'No base surface head model found.'); return;
    end
    HeadModelFile = sStudy.HeadModel(iBase).FileName;

    if mod(size(HeadModel.Gain,2), 3) ~= 0
        bst_report('Error', sProcess, sInputs, ...
            'Base head model must be unconstrained (Gain [nCh x 3*nVert]).'); return;
    end

    % Dirac eigenbasis on the surface (compute if absent / different Tau or K)
    T = in_tess_bst(HeadModel.SurfaceFile, 0);
    needCompute = ~isfield(T,'DiracEigen') || isempty(T.DiracEigen) ...
        || ~isequal([T.DiracEigen.Tau],[Tau Tau]) || any([T.DiracEigen.nModes] ~= Kbasis);
    if needCompute
        bst_progress('text', 'Computing Dirac eigenbasis...');
        DiracEigen = tess_dirac_eigenmodes(HeadModel.SurfaceFile, 'Tau', Tau, 'K', Kbasis);
    else
        DiracEigen = T.DiracEigen;
    end

    nVertHM = size(HeadModel.Gain, 2) / 3;
    nVertDE = sum(arrayfun(@(d) numel(d.GlobalVertices), DiracEigen));
    if nVertHM ~= nVertDE
        bst_report('Error', sProcess, sInputs, sprintf( ...
            'Head model has %d vertices but Dirac eigenbasis covers %d.', nVertHM, nVertDE)); return;
    end

    % Compose (all Kbasis modes) + history
    CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen);
    CompHM = bst_history('add', CompHM, 'dirac_eigenmode_leadfield', ...
        sprintf('Composed Dirac eigenmode leadfield: %d modes (tau=%.3g) from %s', ...
        CompHM.nModes, Tau, HeadModelFile));

    % Save + register (mirrors process_eigenmode_leadfield)
    StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
    OutputFile = bst_process('GetNewFilename', StudyDir, 'headmodel_dirac_eigenmode');
    bst_save(OutputFile, CompHM, 'v7');

    sHeadModel = db_template('headmodel');
    sHeadModel.FileName      = file_short(OutputFile);
    sHeadModel.Comment       = CompHM.Comment;
    sHeadModel.HeadModelType = CompHM.HeadModelType;
    iHM = length(sStudy.HeadModel) + 1;
    sStudy.HeadModel(iHM) = sHeadModel;
    sStudy.iHeadModel     = iHM;
    bst_set('Study', iStudy, sStudy);
    panel_protocols('UpdateNode', 'Study', iStudy);
    panel_protocols('SelectNode', [], file_short(OutputFile));
    db_save();

    % DB side-effect; pass inputs through unchanged.
    OutputFiles = {sInputs.FileName};
end
