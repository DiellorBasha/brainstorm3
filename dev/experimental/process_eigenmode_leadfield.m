function varargout = process_eigenmode_leadfield( varargin )
% PROCESS_EIGENMODE_LEADFIELD: Compose a base head model into the eigenmode basis.
%
% Produces a headmodel_eigenmode_*.mat node whose Gain = L*Phi (the eigenmode
% leadfield). Requires a surface head model and precomputed surface eigenmodes.
%
% Authors: Diellor Basha, 2026
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Compute eigenmode leadfield';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    sProcess.options.nmodes.Comment = 'Number of eigenmodes (0 = all available): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Composes the standard leadfield with the ' ...
        'surface eigenmodes (L*Phi).<BR>Requires a surface head model and precomputed eigenmodes.</FONT>'];
    sProcess.options.label_info.Type = 'label';
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    n = sProcess.options.nmodes.Value{1};
    if n > 0; Comment = sprintf('Eigenmode leadfield (%d modes)', n);
    else;     Comment = 'Eigenmode leadfield (all modes)'; end
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    nModes = sProcess.options.nmodes.Value{1};

    [sStudy, iStudy] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.HeadModel)
        bst_report('Error', sProcess, sInputs, 'No head model available for this study.'); return;
    end

    % Choose a BASE (non-eigenmode) surface head model. The active model may already
    % be an eigenmode leadfield (e.g. when re-running this process); never compose on
    % top of one. Prefer the active head model when it is itself a base surface model.
    iBase = [];
    HeadModel = [];
    candidates = [sStudy.iHeadModel, setdiff(1:numel(sStudy.HeadModel), sStudy.iHeadModel)];
    for ic = candidates
        if ic < 1 || ic > numel(sStudy.HeadModel); continue; end
        try
            hmC = in_bst_headmodel(sStudy.HeadModel(ic).FileName, 0);
        catch
            continue;
        end
        isEig = isfield(hmC,'isEigenmode') && ~isempty(hmC.isEigenmode) && hmC.isEigenmode;
        if ~isEig && strcmpi(hmC.HeadModelType, 'surface')
            iBase = ic; HeadModel = hmC; break;
        end
    end
    if isempty(iBase)
        bst_report('Error', sProcess, sInputs, ...
            'No base surface head model found (only eigenmode leadfields present).'); return;
    end
    HeadModelFile = sStudy.HeadModel(iBase).FileName;

    % Load eigenmodes from the head model's surface
    [Eigenmodes, isComputed] = in_tess_eigenmodes(HeadModel.SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, sInputs, ...
            ['No eigenmodes on surface: ' HeadModel.SurfaceFile '. Run "Compute eigenmodes" first.']); return;
    end

    % Vertex-count consistency (head model is unconstrained: 3 cols/vertex)
    nVertHM = size(HeadModel.Gain, 2) / 3;
    if nVertHM ~= size(Eigenmodes.Vectors, 1)
        bst_report('Error', sProcess, sInputs, sprintf( ...
            ['Head model has %d vertices but eigenmodes have %d.\nRecompute the head model ' ...
             '(computing eigenmodes may have repaired the mesh).'], nVertHM, size(Eigenmodes.Vectors,1))); return;
    end

    % Compose
    CompHM = bst_eigenmode_leadfield(HeadModel, Eigenmodes, 'nModes', nModes);
    CompHM = bst_history('add', CompHM, 'eigenmode_leadfield', ...
        sprintf('Composed eigenmode leadfield: %d modes from %s', CompHM.nModes, HeadModelFile));

    % Save as a new head model node
    StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
    OutputFile = bst_process('GetNewFilename', StudyDir, 'headmodel_eigenmode');
    bst_save(OutputFile, CompHM, 'v7');

    % Register in the DB
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

    % The eigenmode head model is a DB side-effect; pass the input files through
    % unchanged so downstream pipeline steps stay valid.
    OutputFiles = {sInputs.FileName};
end
