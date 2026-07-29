function varargout = process_source_flow_functional( varargin )
% PROCESS_SOURCE_FLOW_FUNCTIONAL: Whole-brain or per-scout flow functionals of an unconstrained
% source (energy |J|^2, enstrophy omega^2, helicity J.(curl J)) as a matrix time series.
%
% USAGE:  OutputFiles = process_source_flow_functional('Run', sProcess, sInputs)
%
% The functionals are BILINEAR in the source J. They are NOT linear maps (so, unlike the flow
% MAPS, they are not a single imaging kernel), but because J = K*b is linear in the sensors they
% ARE a fixed SENSOR-SPACE GRAM form: for a differential operator D and area weights W,
%     f(t) = (D J)^T W (D J) = b(t)^T [ (D K)^T W (D K) ] b(t) = b(t)^T Q b(t),
% with Q an [nChan x nChan] matrix computed ONCE (D K is the fused flow kernel; curl is linear so
% curl(K b) = (curl K) b). The whole-recording series is then b^T Q b per frame -- a tiny 270x270
% quadratic form -- WITHOUT ever reconstructing the 3*nVertices source. This is ~2 orders of
% magnitude cheaper than materializing J every frame and cannot run out of memory (Q is nChan^2).
% Matches the Gram form of the standalone +flow package. A full (ImageGridAmp, no-kernel) source
% has no sensor representation, so it falls back to materialize-and-integrate. Geometry (areas,
% curl) comes from the Covariant operator node.

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
        Rk = in_bst_results(sInputs(iInput).FileName, 0, ...
            'ImagingKernel','GoodChannel','DataFile','nComponents','SurfaceFile','Comment','Time');
        hasKernel = ~isempty(Rk.ImagingKernel);
        if hasKernel
            SurfaceFile = Rk.SurfaceFile;  Comment = Rk.Comment;  nComp = Rk.nComponents;
        else
            R = in_bst_results(sInputs(iInput).FileName, 1);       % full field: materialize (fallback)
            SurfaceFile = R.SurfaceFile;  Comment = R.Comment;  nComp = R.nComponents;
        end
        if nComp ~= 3
            bst_report('Error', sProcess, sInputs(iInput), 'Input is not an unconstrained (3-component) source model.');
            continue;
        end

        Cov = tess_operators(SurfaceFile, 'Covariant');
        w   = VertexAreas(Cov);                                   % [nVtot x 1] lumped vertex areas
        if isempty(AtlasList)
            scoutVerts = {};  Desc = {'whole-brain'};
        else
            [scoutVerts, scoutNames] = GetScoutVertices(AtlasList, SurfaceFile);  Desc = scoutNames(:);
        end

        if hasKernel
            % ---- GRAM FORM: f(t) = b' Q b, no source reconstruction ----
            Qs = BuildGrams(Func, Rk.ImagingKernel, Cov, w, scoutVerts);   % [nCh x nCh] per output row
            sData = in_bst_data(Rk.DataFile, 'F');
            if isstruct(sData.F)                                   % raw continuous -> stream sensors
                ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', Rk.DataFile));
                [Value, Time] = StreamGram(Qs, sData.F, ChannelMat, Rk.GoodChannel, sProcess);
            else                                                   % imported sensors in memory
                Value = ApplyGrams(Qs, sData.F(Rk.GoodChannel, :));  Time = Rk.Time;
            end
        else
            % ---- fallback: full source, no kernel -> materialize + integrate ----
            dens  = Density(Func, R.ImageGridAmp, Cov);
            Value = ReduceDensity(dens, w, scoutVerts);            Time = R.Time;
        end

        % Save as a matrix time series
        FileMat = db_template('matrixmat');
        FileMat.Value        = Value;
        FileMat.Time         = Time;
        FileMat.Description   = Desc;
        FileMat.Comment      = [Comment ' | flow:' Func];
        FileMat.DisplayUnits  = 'flow';
        FileMat = bst_history('add', FileMat, 'flow', ['Flow functional (Gram): ' Func]);
        sStudy = bst_get('Study', sInputs(iInput).iStudy);
        OutputFiles{iInput} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['matrix_flow_' Func]);
        bst_save(OutputFiles{iInput}, FileMat, 'v6');
        db_add_data(sInputs(iInput).iStudy, OutputFiles{iInput}, FileMat);
    end
end


%% ===== BUILD SENSOR-SPACE GRAM(S)  Q = sym( (D K)' * (W .* (E K)) )  [nCh x nCh] =====
% The functional f = (D J)' W (E J) with J = K b becomes f = b' Q b, Q = sym( AK' (W .* BK) ).
% AK/BK are the operator-applied kernels (computed ONCE); per-scout Grams reuse them with a masked
% weight, so curl(K) is never recomputed. energy: A=B=K; enstrophy: A=B=curl K; helicity: A=K, B=curlVec K.
function Qs = BuildGrams(Func, K, Cov, w, scoutVerts)
    switch Func
        case 'energy'
            AK = K;                                   BK = K;                    scalarWeight = false;
        case 'enstrophy'
            AK = bst_curl(K, [], 'Ambient', [], Cov); BK = AK;                   scalarWeight = true;   % [nV x nCh]
        case 'helicity'
            AK = K;                                   BK = CurlVector(K, Cov);   scalarWeight = false;  % [3nV x nCh]
        otherwise
            error('Unknown functional "%s".', Func);
    end
    if isempty(scoutVerts), weightSets = {w};
    else
        weightSets = cell(numel(scoutVerts), 1);
        for k = 1:numel(scoutVerts)
            ws = zeros(numel(w), 1);  ws(scoutVerts{k}) = w(scoutVerts{k});  weightSets{k} = ws;
        end
    end
    Qs = cell(numel(weightSets), 1);
    for k = 1:numel(weightSets)
        if scalarWeight, Wt = weightSets{k};                            % weight on per-vertex scalar (nV)
        else,            Wt = reshape(repmat(weightSets{k}(:)', 3, 1), [], 1);  end  % on 3-vector rows (3nV)
        Qraw = AK' * (Wt .* BK);
        Qs{k} = 0.5 * (Qraw + Qraw');                                   % symmetric quadratic form
    end
end


%% ===== STREAM SENSORS AND APPLY THE GRAM(S)  f(t) = b(t)' Q b(t) =====
% Reads the raw sensors (nChan rows only -- tiny) via the canonical reader
% panel_record('ReadRawBlock', ..., 'all', 1), so blocks can be large; no source is materialized.
function [Value, Time] = StreamGram(Qs, sFile, ChannelMat, GoodChannel, sProcess)
    sfreq = sFile.prop.sfreq;  tStart = sFile.prop.times(1);  tStop = sFile.prop.times(2);
    blockSec = 30;                                          % b is nChan rows -> big blocks are cheap
    Value = zeros(numel(Qs), 0);  Time = zeros(1, 0);
    tBlock = tStart;  iBlock = 0;  nBlocks = ceil((tStop - tStart) / blockSec);
    while tBlock < tStop - 0.5/sfreq
        tBlockEnd = min(tBlock + blockSec, tStop);
        [b, T] = panel_record('ReadRawBlock', sFile, ChannelMat, 1, [tBlock tBlockEnd], 0, 1, 'all', 1);
        Value = [Value, ApplyGrams(Qs, b(GoodChannel, :))];  %#ok<AGROW>
        Time  = [Time,  T];                                  %#ok<AGROW>
        tBlock = tBlockEnd + 1/sfreq;  iBlock = iBlock + 1;
        if nargin >= 5 && ~isempty(sProcess)
            bst_progress('text', sprintf('Flow functional (Gram): block %d/%d', iBlock, nBlocks));
        end
    end
end

function Value = ApplyGrams(Qs, b)
    Value = zeros(numel(Qs), size(b, 2));
    for r = 1:numel(Qs)
        Value(r,:) = sum(b .* (Qs{r} * b), 1);               % b' Q b per frame
    end
end


%% ===== FUNCTIONAL DENSITY (per-vertex; fallback path for full sources) =====
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


%% ===== INTEGRATE A PER-VERTEX DENSITY OVER CORTEX OR PER SCOUT (fallback) =====
function Value = ReduceDensity(dens, w, scoutVerts)
    if isempty(scoutVerts)
        Value = sum(w .* dens, 1);
    else
        Value = zeros(numel(scoutVerts), size(dens, 2));
        for k = 1:numel(scoutVerts)
            v = scoutVerts{k};  Value(k,:) = sum(w(v) .* dens(v,:), 1);
        end
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
