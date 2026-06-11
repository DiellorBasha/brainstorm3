function ResultsFile = make_alpha_dirac_source()
% MAKE_ALPHA_DIRAC_SOURCE  Create an unconstrained (nComponents=3) Dirac source
% node for the alpha-band data block (data_block001_band, 7-13 Hz) so the
% quiver/vector-field overlay in figure_3d has a real vector field to draw.
%
% The function runs bst_inverse_dirac on the UNCONSTRAINED surface head model
% (headmodel_surf_os_meg) — bst_inverse_dirac internally computes the Dirac
% eigenbasis transform — and writes the resulting ImagingKernel (3*nVert x nMEG)
% as a new results_*.mat under the same study directory, then registers it in
% the Brainstorm DB.
%
% Running this function a second time creates a second (uniquely-named) node —
% that is intentional; it never deletes existing nodes.
%
% Returns the SHORT (protocol-relative) path of the newly created file.
%
% Authors: Diellor Basha, 2026

    % ------------------------------------------------------------------ %
    % 1.  Locate the study, the alpha-band data, head model, noise cov
    % ------------------------------------------------------------------ %
    cond = 'Subject01/S01_AEF_20131218_01_notch';
    [sStudy, iStudy] = bst_get('StudyWithCondition', cond);
    if isempty(sStudy)
        error('make_alpha_dirac_source:noStudy', ...
            'Condition not found in the current protocol: %s', cond);
    end

    % Alpha-band data node
    iData = find(~cellfun(@isempty, ...
        regexp({sStudy.Data.FileName}, 'data_block001_band', 'once')), 1);
    if isempty(iData)
        error('make_alpha_dirac_source:noAlphaData', ...
            'alpha band data node (data_block001_band) not found in %s', cond);
    end
    DataFile = sStudy.Data(iData).FileName;   % short (protocol-relative)

    % Unconstrained surface head model (headmodel_surf_os_meg, Gain [340 x 61452])
    iHM_surf = find(~cellfun(@isempty, ...
        regexp({sStudy.HeadModel.FileName}, 'headmodel_surf_os_meg', 'once')), 1);
    if isempty(iHM_surf)
        error('make_alpha_dirac_source:noSurfHM', ...
            'headmodel_surf_os_meg not found in %s', cond);
    end
    HMfile = sStudy.HeadModel(iHM_surf).FileName;
    HMos   = in_bst_headmodel(HMfile, 0);   % loads Gain [nCh x 3*nVert]

    % Noise covariance
    NCfile = sStudy.NoiseCov(1).FileName;
    NC     = load(file_fullpath(NCfile));    % has .NoiseCov [340 x 340]

    % Channel info
    Chan = in_bst_channel(sStudy.Channel.FileName);

    % ------------------------------------------------------------------ %
    % 2.  Select MEG channels (same convention as bst_inverse_linear_2018)
    % ------------------------------------------------------------------ %
    G       = double(HMos.Gain);               % [340 x 61452]
    isMEG   = strcmpi({Chan.Channel.Type}', 'MEG');
    iMEG    = find(isMEG(:)');                % [1 x nMEG] column indices into 340

    % Restrict head model and noise cov to MEG rows
    HMf          = HMos;
    HMf.Gain     = G(iMEG, :);               % [nMEG x 61452]

    OPT = struct();
    OPT.NoiseCovMat.NoiseCov = NC.NoiseCov(iMEG, iMEG);  % [nMEG x nMEG]
    OPT.ChannelTypes         = {Chan.Channel(iMEG).Type};
    OPT.NoiseMethod          = 'reg';
    OPT.NoiseReg             = 0.1;
    OPT.SnrMethod            = 'fixed';
    OPT.SnrFixed             = 3;
    OPT.InverseMeasure       = 'amplitude';
    OPT.nModes               = 400;
    OPT.Tau                  = 0.5;

    % ------------------------------------------------------------------ %
    % 3.  Run the Dirac inverse  -> ImagingKernel [3*nVert x nMEG]
    % ------------------------------------------------------------------ %
    fprintf('make_alpha_dirac_source > Running bst_inverse_dirac ...\n');
    R = bst_inverse_dirac(HMf, OPT);
    % R.ImagingKernel is [3*nVert x nMEG] on raw data (whitener already folded in)
    fprintf('make_alpha_dirac_source > ImagingKernel: [%d x %d]\n', ...
        size(R.ImagingKernel, 1), size(R.ImagingKernel, 2));

    % ------------------------------------------------------------------ %
    % 4.  Build the Results struct, using the existing Dirac kernel as a
    %     field-set template (provides ChannelFlag, nAvg, Leff, etc.)
    % ------------------------------------------------------------------ %
    % Find the existing Dirac kernel file to use as struct template
    iDiracKernel = find(~cellfun(@isempty, ...
        regexp({sStudy.Result.FileName}, 'DiracEig_KERNEL', 'once')), 1);
    if isempty(iDiracKernel)
        % Fall back to any non-link results file
        iDiracKernel = find(~cellfun(@(f) startsWith(f,'link|'), ...
            {sStudy.Result.FileName}), 1);
    end
    if isempty(iDiracKernel)
        error('make_alpha_dirac_source:noTemplate', ...
            'No existing results_*.mat to template from in %s', cond);
    end
    Tmpl = load(file_fullpath(sStudy.Result(iDiracKernel).FileName));

    % Load the alpha-band data to obtain its Time vector (needed for display)
    Dmat = in_bst_data(DataFile);             % loads .Time, .ChannelFlag, etc.

    % -- Override all kernel-specific fields -- %
    Tmpl.ImagingKernel = R.ImagingKernel;     % [3*nVert x nMEG]
    Tmpl.ImageGridAmp  = [];
    Tmpl.Std           = [];
    Tmpl.nComponents   = 3;                   % UNCONSTRAINED: 3 components per vertex
    Tmpl.Function      = 'dirac';
    Tmpl.Comment       = 'Dirac vector (amplitude) | alpha band';

    % Timing: copy from the data block so the viewer can animate correctly
    Tmpl.Time          = Dmat.Time;           % [1 x 1501]

    % DataFile links this result to the alpha-band data node
    Tmpl.DataFile      = DataFile;

    % Head model provenance
    Tmpl.HeadModelFile = HMfile;
    Tmpl.HeadModelType = HMos.HeadModelType;  % 'surface'
    Tmpl.SurfaceFile   = HMos.SurfaceFile;    % 'Subject01/tess_cortex_pial_low.mat'

    % Grid geometry (empty for surface kernel, same as template / MNE convention)
    Tmpl.GridLoc       = [];
    Tmpl.GridOrient    = [];
    Tmpl.GridAtlas     = [];

    % ChannelFlag: take from the data file (full 340-channel flag)
    Tmpl.ChannelFlag   = Dmat.ChannelFlag;    % [340 x 1]

    % GoodChannel: indices of the MEG channels used (into the full 340-ch array)
    % These must correspond 1-to-1 with columns of ImagingKernel.
    Tmpl.GoodChannel   = iMEG;               % [1 x nMEG] e.g. [1 x 274]

    % Inherit nAvg and Leff from the data node
    Tmpl.nAvg          = Dmat.nAvg;
    Tmpl.Leff          = Dmat.Leff;

    % Carry through Dirac-specific metadata that bst_inverse_dirac returned
    diracFields = {'Eigenvalues','ModeHemisphere','DiracEigenFile','DiracTau','ImagingKernelMode'};
    for k = 1:numel(diracFields)
        f = diracFields{k};
        if isfield(R, f)
            Tmpl.(f) = R.(f);
        end
    end

    % ------------------------------------------------------------------ %
    % 5.  Write to a new uniquely-timestamped file and register in DB
    % ------------------------------------------------------------------ %
    FullDataFile = file_fullpath(DataFile);
    DataDir      = fileparts(FullDataFile);
    ResultsFile  = bst_process('GetNewFilename', DataDir, 'results_dirac_vec');

    bst_save(ResultsFile, Tmpl, 'v6');
    db_add_data(iStudy, ResultsFile, Tmpl);
    db_reload_studies(iStudy);

    ResultsFile = file_short(ResultsFile);
    fprintf('make_alpha_dirac_source > Created unconstrained Dirac source node:\n  %s\n', ResultsFile);
end
