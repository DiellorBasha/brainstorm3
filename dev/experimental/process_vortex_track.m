function varargout = process_vortex_track( varargin )
% PROCESS_VORTEX_TRACK  Track Dirac vortex cores (and sources/sinks) over time.
% Detects persistence-ranked cores per frame (bst_helmholtz) and links them
% into trajectories (bst_vortex_track), saved as a dipoles file (view_dipoles).
% @=============================================================================
% Author: Diellor Basha, 2026
    eval(macro_method);
end

%% ===== DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Track vortex cores (Dirac)';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 0;
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'dipoles'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    sProcess.options.minpers.Comment = 'Min persistence (0 = keep all): ';
    sProcess.options.minpers.Type    = 'value';
    sProcess.options.minpers.Value   = {0, '', 6};
    sProcess.options.maxjump.Comment = 'Max core jump per frame: ';
    sProcess.options.maxjump.Type    = 'value';
    sProcess.options.maxjump.Value   = {10, 'mm', 1};
    sProcess.options.tracksrc.Comment = 'Also track sources/sinks (phi)';
    sProcess.options.tracksrc.Type    = 'checkbox';
    sProcess.options.tracksrc.Value   = 1;
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {};
    MinPers = sProcess.options.minpers.Value{1};
    MaxJump = sProcess.options.maxjump.Value{1} / 1000;     % mm -> m
    TrackSrc = sProcess.options.tracksrc.Value;
    TimeWindow = sProcess.options.timewindow.Value{1};

    % --- source file (must be unconstrained / 3-component) ---
    sRes = in_bst_results(sInput.FileName, 0);
    if isempty(sRes.nComponents) || (sRes.nComponents ~= 3)
        bst_report('Error', sProcess, sInput, 'Vortex tracking requires an unconstrained (3-component) Dirac source.');
        return;
    end
    % data over the window
    if isempty(sRes.DataFile), DataMat.Time = sRes.Time; DataMat.F = [];
    else, DataMat = in_bst_data(sRes.DataFile); end
    if isempty(sRes.ImageGridAmp)
        F = DataMat.F(sRes.GoodChannel, :);
        J = sRes.ImagingKernel * F;
    else
        J = sRes.ImageGridAmp;
    end
    Time = sRes.Time; if isempty(Time) || numel(Time)~=size(J,2), Time = DataMat.Time; end
    if ~isempty(TimeWindow)
        sb = bst_closest(TimeWindow, Time); iW = sb(1):sb(2);
    else
        iW = 1:size(J,2);
    end
    J = J(:, iW); tW = Time(iW);

    % --- operators + per-frame decomposition ---
    SurfaceFile = sRes.SurfaceFile;
    Surf = in_tess_bst(SurfaceFile, 0);
    Mani = tess_manifold(SurfaceFile);
    Cov = tess_operators(SurfaceFile, 'Covariant');  LBO = tess_operators(SurfaceFile, 'Laplace-Beltrami');
    Op = bst_helmholtz('Prepare', {Cov, LBO}, Mani, Surf, 'Domain','vertex');
    nT = numel(tW);
    coresV = cell(1,nT); coresS = cell(1,nT);
    bst_progress('start', 'Vortex tracking', 'Decomposing frames...', 1, nT);
    for t = 1:nT
        Ht = bst_helmholtz('Frame', Op, J(:,t));
        coresV{t} = Ht.Cores;  coresS{t} = Ht.Sources;
        bst_progress('inc', 1);
    end
    bst_progress('stop');

    % --- link ---
    Tv = bst_vortex_track(coresV, 'MinPersistence', MinPers, 'MaxJump', MaxJump);
    Ts = i_empty_tracks();
    if TrackSrc, Ts = bst_vortex_track(coresS, 'MinPersistence', MinPers, 'MaxJump', MaxJump); end
    if isempty(Tv) && isempty(Ts)
        bst_report('Warning', sProcess, sInput, 'No vortex cores found in the selected window.');
    end

    % --- assemble dipoles ---
    VN = Surf.VertNormals;
    DipolesMat = i_tracks_to_dipoles(Tv, Ts, tW, VN);
    DipolesMat.Comment  = sprintf('Vortex tracks (%d+%d)', numel(Tv), numel(Ts));
    DipolesMat.DataFile = sInput.FileName;
    DipolesMat = bst_history('add', DipolesMat, 'vortextrack', ['Generated from: ' sInput.FileName]);

    % --- save + db ---
    [sStudy, iStudy] = bst_get('AnyFile', sInput.FileName);
    ProtocolInfo = bst_get('ProtocolInfo');
    DipoleFile = file_unique(bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sStudy.FileName), 'dipoles_vortextrack.mat'));
    bst_save(DipoleFile, DipolesMat);
    Bst = db_template('Dipoles');
    Bst.FileName = file_short(DipoleFile);
    Bst.Comment  = DipolesMat.Comment;
    Bst.DataFile = sInput.FileName;
    sStudy.Dipoles(end+1) = Bst;
    bst_set('Study', iStudy, sStudy);
    panel_protocols('UpdateNode', 'Study', iStudy);
    db_save();
    OutputFiles{1} = DipoleFile;
end

%% ===== helpers =====
function s = i_empty_tracks()
    s = struct('frames',{},'iVertex',{},'pos',{},'persistence',{}, ...
               'chirality',{},'birthFrame',{},'deathFrame',{});
end

function D = i_tracks_to_dipoles(Tv, Ts, tW, VN)
% Each track -> one dipole .Index group; Loc=pos, Amplitude=chirality*normPers*normal.
    fp = [];
    for k = 1:numel(Tv), p = Tv(k).persistence; fp = [fp, p(isfinite(p))]; end %#ok<AGROW>
    for k = 1:numel(Ts), p = Ts(k).persistence; fp = [fp, p(isfinite(p))]; end %#ok<AGROW>
    maxP = max([fp, eps]);
    D = struct(); D.Time = tW(:)'; D.Dipole = repmat(i_empty_dip(), 0, 1);
    names = {};
    gi = 0;
    [D, gi, names] = i_add_set(D, gi, names, Tv, tW, VN, maxP, 'Vortex');
    [D, gi, names] = i_add_set(D, gi, names, Ts, tW, VN, maxP, 'Source');  %#ok<ASGLU>
    D.DipoleNames = names;
    if isempty(D.Dipole), D.Subset = []; else, D.Subset = unique([D.Dipole.Index]); end
end

function d = i_empty_dip()
    d = struct('Index',0,'Time',0,'Origin',[0 0 0],'Loc',[0;0;0],'Amplitude',[0;0;0], ...
               'Goodness',0,'Errors',0,'Noise',[],'SingleError',[],'ErrorMatrix',[], ...
               'ConfVol',[],'Probability',[],'NoiseEstimate',[],'Perform',0);
end

function [D, gi, names] = i_add_set(D, gi, names, T, tW, VN, maxP, kind)
    for k = 1:numel(T)
        gi = gi + 1;  tr = T(k);
        if strcmp(kind,'Vortex'); lbl = sprintf('Vortex%s #%d', i_sgn(tr.chirality), gi);
        else;                     lbl = sprintf('%s #%d', i_srcname(tr.chirality), gi); end
        names{gi} = lbl; %#ok<AGROW>
        for j = 1:numel(tr.frames)
            d = i_empty_dip();
            d.Index = gi;  d.Time = tW(tr.frames(j));
            d.Loc = tr.pos(j,:)';
            pn = tr.persistence(j); if ~isfinite(pn), pn = maxP; end
            d.Amplitude = (tr.chirality * (pn/maxP) * VN(tr.iVertex(j),:))';
            d.Goodness = min(pn, maxP)/maxP;  d.Perform = pn;
            D.Dipole(end+1) = d; %#ok<AGROW>
        end
    end
end

function s = i_sgn(ch),     if ch>=0, s='+'; else, s='-'; end, end
function s = i_srcname(ch), if ch>=0, s='Source'; else, s='Sink'; end, end
