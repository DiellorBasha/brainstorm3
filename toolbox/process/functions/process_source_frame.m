function varargout = process_source_frame( varargin )
% PROCESS_SOURCE_FRAME: Opt-in whole-series frame scalogram (source).
%
% Batch (never interactive) version of the dynamics differential-frame overlay.
% Computes an itersine tight-frame scalogram of a source result over the WHOLE
% series, in contiguous non-overlapping windows (never loading the whole series
% into memory at once -- each window is paged in via bst_memory), and saves the
% concatenated result as a single whole-series scalogram (Timefreq file, 3 rows
% {Global,LH,RH}, full time axis).
%
% Reuses Task 1's bst_eigenwavelet('Scalogram', ax, gCell, C) and Task 1's
% itersine kernel (bst_eigfilter_design_itersine); the frame + per-window
% projection are built inline (mirrors bst_dynamics/panel_bst_dynamics, but
% those are the interactive single-frame GUI path -- this process is the
% opt-in whole-series batch path).
%
% Input  : a source result file (results).
% Output : a whole-series scalogram (timefreq) file.
%
% USAGE:  OutputFiles = process_source_frame('Run', sProcess, sInputs)
%
% Authors: Diellor Basha, 2026

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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    sProcess.Comment     = 'Frame scalogram (source)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 339;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm';
    % Input/Output
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === NUMBER OF FRAME MEMBERS ===
    sProcess.options.nframe.Comment = 'Itersine frame members: ';
    sProcess.options.nframe.Type    = 'value';
    sProcess.options.nframe.Value   = {6, '', 0};

    % === EIGENBASIS VARIANT ===
    sProcess.options.variant.Comment = 'Eigenbasis variant: ';
    sProcess.options.variant.Type    = 'combobox_label';
    sProcess.options.variant.Value   = {'Geometric', {'Geometric', 'Connectomic'; 'Geometric', 'Connectomic'}};

    % === WINDOW LENGTH ===
    sProcess.options.winsec.Comment = 'Window length (s): ';
    sProcess.options.winsec.Type    = 'value';
    sProcess.options.winsec.Value   = {4, 's', 2};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sprintf('Frame scalogram (source): %d members', sProcess.options.nframe.Value{1});
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    N   = sProcess.options.nframe.Value{1};
    varOpt = sProcess.options.variant.Value;
    variant = 'Laplace-Beltrami'; if iscell(varOpt) && ~isempty(varOpt) && strcmpi(varOpt{1},'Connectomic'), variant = 'LB-Connectome'; end
    winsec = sProcess.options.winsec.Value{1};
    for iIn = 1:numel(sInputs)
        R = in_bst_results(sInputs(iIn).FileName, 0, 'SurfaceFile','Time','DataFile');
        ax = bst_eigen('Axes', struct('SurfaceFile',R.SurfaceFile, 'Variant',variant, 'nModes',60, ...
                        'TimeWindow',[R.Time(1) R.Time(min(2,end))], 'SampleRate', 1/median(diff(R.Time))));
        lmax = max(ax.Lambda{1}(:));
        gCell = cell(1,N); for ii=1:N, gCell{ii} = bst_eigfilter_design_itersine(struct('member',ii,'Nf',N,'lmax',lmax)); end
        nV = 0; for h=1:numel(ax.GlobalVertices), nV=max(nV,max(ax.GlobalVertices{h}(:))); end
        tv = R.Time;  Fs = 1/median(diff(tv));  nWin = max(1, round(winsec*Fs));
        starts = 1:nWin:numel(tv);
        [iDS, iRes] = bst_memory('LoadResultsFileFull', sInputs(iIn).FileName);   % load once; page windows below
        E = [];  Res = [];  Tall = [];
        for s0 = starts
            iWin = s0:min(s0+nWin-1, numel(tv));
            F = double(bst_memory('GetResultsValues', iDS, iRes, [], iWin, 0));   % reconstruct only this window
            % reduce to scalar magnitude, project per hemi
            F = i_reduce(F, nV);
            C = cell(1,numel(ax.Phi));
            for h=1:numel(ax.Phi), gv=ax.GlobalVertices{h}(:); C{h}=manifold_ft(ax.Phi{h}, ax.Mass{h}, F(gv,:)); end
            scal = bst_eigenwavelet('Scalogram', ax, gCell, C);
            E = cat(2, E, scal.energy);  Res = [Res, scal.residual];  Tall = [Tall, tv(iWin)]; %#ok<AGROW>
        end
        FileMat = db_template('timefreqmat');
        FileMat.TF = E;  FileMat.Time = Tall;  FileMat.Freqs = scal.centers(:);
        FileMat.RowNames = {'Global','LH','RH'};  FileMat.Measure='power';  FileMat.Method='framescalogram';
        FileMat.DataType='matrix';  FileMat.SurfaceFile=R.SurfaceFile;  FileMat.nAvg=1; FileMat.Leff=1;
        FileMat.DataFile = file_short(sInputs(iIn).FileName);
        FileMat.Comment = sprintf('Frame scalogram (series) | itersine x%d, %s', N, variant);
        OutFile = bst_process('GetNewFilename', bst_fileparts(file_fullpath(sInputs(iIn).FileName)), 'timefreq_framescalo_series');
        bst_save(OutFile, FileMat, 'v6');  db_add_data(sInputs(iIn).iStudy, file_short(OutFile), FileMat);
        OutputFiles{end+1} = OutFile; %#ok<AGROW>
    end
end

% scalar magnitude reduction (k*nV -> nV), mirrors the panel's i_paintable_scalar
function s = i_reduce(F, nV)
    if ~isreal(F), F = abs(F); end
    if size(F,1)==nV, s = F; return; end
    if mod(size(F,1),nV)==0, nc=size(F,1)/nV; s = reshape(sqrt(sum(reshape(F,nc,nV,[]).^2,1)),nV,[]); else, s = F; end
end
