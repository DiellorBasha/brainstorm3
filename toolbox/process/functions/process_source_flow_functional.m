function varargout = process_source_flow_functional( varargin )
% PROCESS_SOURCE_FLOW_FUNCTIONAL: Whole-brain or per-scout flow functionals of an unconstrained
% source (energy |J|^2, enstrophy omega^2, helicity J.(curl J)) as a matrix time series.
%
% USAGE:  OutputFiles = process_source_flow_functional('Run', sProcess, sInputs)
%
% The functionals are QUADRATIC in the source, so unlike the flow MAPS (process_source_flow)
% they cannot be kernels: the source is materialized and the functional integrated over the
% cortex (or per scout), giving a scalar / per-scout time series (a matrix file) that PSD /
% Welch analyze. Geometry (areas, curl) comes from the Covariant operator node.

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
    sProcess.Comment     = 'Flow functional (energy/enstrophy/helicity)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 339;
    sProcess.Description = '';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.options.label1.Comment = ['Integrated flow functional of an unconstrained source,<BR>' 10 ...
        'as a time series for PSD / Welch. Whole-brain, or one row<BR>' 10 ...
        'per scout if an atlas is selected.<BR><BR>Functional:'];
    sProcess.options.label1.Type = 'label';
    sProcess.options.func.Comment = {'<B>Energy</B> |J|^2', '<B>Enstrophy</B> omega^2', '<B>Helicity</B> J.(curl J)'; ...
        'energy', 'enstrophy', 'helicity'};
    sProcess.options.func.Type  = 'radio_label';
    sProcess.options.func.Value = 'energy';
    sProcess.options.scouts.Comment = 'Per-scout (leave empty = whole-brain):';
    sProcess.options.scouts.Type    = 'scout';
    sProcess.options.scouts.Value    = {};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    Func = sProcess.options.func.Value;
    AtlasList = sProcess.options.scouts.Value;
    for iInput = 1:numel(sInputs)
        % Resolve the source WITHOUT materializing: a kernel over a continuous (raw) file is
        % streamed block-by-block; anything else (short imported source) is materialized.
        Rk = in_bst_results(sInputs(iInput).FileName, 0, ...
            'ImagingKernel','GoodChannel','DataFile','nComponents','SurfaceFile','Comment','Time');
        isRawKernel = ~isempty(Rk.ImagingKernel) && ~isempty(Rk.DataFile);
        if isRawKernel
            sRawData = in_bst_data(Rk.DataFile, 'F');
            isRawKernel = isstruct(sRawData.F);                 % a continuous-file descriptor => raw link
        end

        if isRawKernel
            SurfaceFile = Rk.SurfaceFile;  Comment = Rk.Comment;  nComp = Rk.nComponents;
        else
            R = in_bst_results(sInputs(iInput).FileName, 1);    % materialize (functional is quadratic)
            SurfaceFile = R.SurfaceFile;  Comment = R.Comment;  nComp = R.nComponents;
        end
        if nComp ~= 3
            bst_report('Error', sProcess, sInputs(iInput), 'Input is not an unconstrained (3-component) source model.');
            continue;
        end

        Cov = tess_operators(SurfaceFile, 'Covariant');
        w   = VertexAreas(Cov);                                 % [nVtot x 1] lumped vertex areas
        if isempty(AtlasList)
            scoutVerts = {};  Desc = {'whole-brain'};
        else
            [scoutVerts, scoutNames] = GetScoutVertices(AtlasList, SurfaceFile);  Desc = scoutNames(:);
        end

        if isRawKernel
            [Value, Time] = StreamFunctional(Func, Rk, sRawData.F, Cov, w, scoutVerts, sProcess);
        else
            dens  = Density(Func, R.ImageGridAmp, Cov);         % [nVtot x nT] per-vertex density
            Value = ReduceDensity(dens, w, scoutVerts);
            Time  = R.Time;
        end

        % Save as a matrix time series
        FileMat = db_template('matrixmat');
        FileMat.Value        = Value;
        FileMat.Time         = Time;
        FileMat.Description   = Desc;
        FileMat.Comment      = [Comment ' | flow:' Func];
        FileMat.DisplayUnits  = 'flow';
        FileMat = bst_history('add', FileMat, 'flow', ['Flow functional: ' Func]);
        sStudy = bst_get('Study', sInputs(iInput).iStudy);
        OutputFiles{iInput} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['matrix_flow_' Func]);
        bst_save(OutputFiles{iInput}, FileMat, 'v6');
        db_add_data(sInputs(iInput).iStudy, OutputFiles{iInput}, FileMat);
    end
end


%% ===== STREAM A FUNCTIONAL OVER A RAW KERNEL-LINK (no full-source materialization) =====
% Reads the raw in memory-bounded time blocks via the canonical reader
% (panel_record('ReadRawBlock', ..., RemoveBaseline='all', UseSsp=1) -- VERIFIED bit-exact against
% bst_process('LoadInputFile')), applies the kernel per block, integrates the functional over the
% cortex (or per scout), and concatenates the scalar/per-scout time series.
function [Value, Time] = StreamFunctional(Func, Rk, sFile, Cov, w, scoutVerts, sProcess)
    K  = Rk.ImagingKernel;  GoodChannel = Rk.GoodChannel;
    ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', Rk.DataFile));
    sfreq = sFile.prop.sfreq;  tStart = sFile.prop.times(1);  tStop = sFile.prop.times(2);
    % memory-bounded block length: cap the per-block source array (~1 GB) since it is 3*nV x nT
    nSources   = size(K, 1);
    maxSamples = max(round(2*sfreq), floor(1e9 / (nSources * 8)));   % >=2 s worth of columns budget
    blockSec   = max(1, maxSamples / sfreq);
    nRows = max(1, numel(scoutVerts));
    Value = zeros(nRows, 0);  Time = zeros(1, 0);
    tBlock = tStart;  iBlock = 0;  nBlocks = ceil((tStop - tStart) / blockSec);
    while tBlock < tStop - 0.5/sfreq
        tBlockEnd = min(tBlock + blockSec, tStop);
        [F, T] = panel_record('ReadRawBlock', sFile, ChannelMat, 1, [tBlock tBlockEnd], 0, 1, 'all', 1);
        Jblock = K * F(GoodChannel, :);                          % [3nV x nT_block]
        dens   = Density(Func, Jblock, Cov);                     % [nV x nT_block]
        Value  = [Value, ReduceDensity(dens, w, scoutVerts)];    %#ok<AGROW>
        Time   = [Time,  T];                                     %#ok<AGROW>
        tBlock = tBlockEnd + 1/sfreq;
        iBlock = iBlock + 1;
        if nargin >= 7 && ~isempty(sProcess)
            bst_progress('text', sprintf('Flow functional: block %d/%d', iBlock, nBlocks));
        end
    end
end


%% ===== INTEGRATE A PER-VERTEX DENSITY OVER CORTEX OR PER SCOUT =====
function Value = ReduceDensity(dens, w, scoutVerts)
    if isempty(scoutVerts)
        Value = sum(w .* dens, 1);                               % [1 x nT] whole-brain
    else
        Value = zeros(numel(scoutVerts), size(dens, 2));
        for k = 1:numel(scoutVerts)
            v = scoutVerts{k};  Value(k,:) = sum(w(v) .* dens(v,:), 1);
        end
    end
end


%% ===== FUNCTIONAL DENSITY (per-vertex; caller applies area weight) =====
function d = Density(Func, J, Cov)
    Jx = J(1:3:end,:);  Jy = J(2:3:end,:);  Jz = J(3:3:end,:);
    switch Func
        case 'energy'
            d = Jx.^2 + Jy.^2 + Jz.^2;
        case 'enstrophy'
            om = bst_curl(J, [], 'Ambient', [], Cov);           % [nVtot x nT] scalar vorticity
            d = om.^2;
        case 'helicity'
            cv = CurlVector(J, Cov);                             % [3nVtot x nT]
            d = Jx.*cv(1:3:end,:) + Jy.*cv(2:3:end,:) + Jz.*cv(3:3:end,:);
        otherwise
            error('Unknown functional "%s".', Func);
    end
end


%% ===== AMBIENT CURL VECTOR (grad x J), area-weighted to vertices =====
function cv = CurlVector(J, Cov)
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    cv = zeros(3*nVtot, size(J,2));
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));  nFh = size(C.Faces,1);
        Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
        Wfv = bst_face2vertex(C.Faces, C.FaceArea);
        Jx = J(3*(vH-1)+1,:);  Jy = J(3*(vH-1)+2,:);  Jz = J(3*(vH-1)+3,:);   % [nVh x nT] vertex field
        cvx = Gy*Jz - Gz*Jy;  cvy = Gz*Jx - Gx*Jz;  cvz = Gx*Jy - Gy*Jx;      % (grad x J) per face [nFh x nT]
        cv(3*(vH-1)+1,:) = Wfv*cvx;  cv(3*(vH-1)+2,:) = Wfv*cvy;  cv(3*(vH-1)+3,:) = Wfv*cvz;
    end
end


%% ===== LUMPED VERTEX AREAS [nVtot x 1] =====
function w = VertexAreas(Cov)
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    w = zeros(nVtot,1);
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));
        a = accumarray(C.Faces(:), repmat(C.FaceArea,3,1), [numel(vH) 1]) / 3;
        w(vH) = a;
    end
end


%% ===== DIRECT ENERGY SERIES (for validation) =====
function E = EnergySeries(J, w)
    Jx = J(1:3:end,:);  Jy = J(2:3:end,:);  Jz = J(3:3:end,:);
    E = sum(w .* (Jx.^2 + Jy.^2 + Jz.^2), 1);
end


%% ===== SCOUT VERTEX LISTS (from an atlas-list option {AtlasName,{Labels}}) =====
function [verts, names] = GetScoutVertices(AtlasList, SurfaceFile)
    sSurf = in_tess_bst(SurfaceFile);
    AtlasName = AtlasList{1};  Labels = AtlasList{2};
    ia = find(strcmp({sSurf.Atlas.Name}, AtlasName), 1);
    if isempty(ia), error('Atlas "%s" not found on the surface.', AtlasName); end
    scouts = sSurf.Atlas(ia).Scouts;
    verts = {};  names = {};
    for k = 1:numel(Labels)
        is = find(strcmp({scouts.Label}, Labels{k}), 1);
        if isempty(is), continue; end
        verts{end+1} = double(scouts(is).Vertices(:));  %#ok<AGROW>
        names{end+1} = scouts(is).Label;                %#ok<AGROW>
    end
end
