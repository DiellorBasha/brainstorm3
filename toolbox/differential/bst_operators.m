function [OutputFiles, Messages, isError] = bst_operators(Data, OPTIONS)
% BST_OPERATORS: Differential-geometry analysis of surface-mapped data on the cortical MANIFOLD
%   domain (the operator_ / manifold_ files). The manifold-domain sibling of BST_EIGEN.
%
% USAGE:  [OutputFiles, Messages, isError] = bst_operators(Data, OPTIONS)
%                                   OPTIONS = bst_operators();
%
% DESCRIPTION:
%     Where bst_eigen analyzes a source map in the SPECTRAL eigenbasis (eigen_ file, the
%     spatial-frequency axis), bst_operators applies discrete DIFFERENTIAL operators read from
%     the operator_ / manifold_ nodes (nxr-built, cached on disk) directly on the cortical
%     manifold:
%
%       bst_eigen     : source map + eigenbasis (Lambda, Phi)  -> spectral transform
%       bst_operators : source map + operators (grad, L, M, D) -> differential transform
%
%     It has NO temporal sibling (cf. bst_eigen <-> bst_timefreq): time-resolved differential
%     analysis is driven by process_ wrappers that loop this orchestrator. I/O-free engines:
%     the nodes are resolved/loaded here, the per-Method engines never touch the DB or disk.
%
%     GEOMETRY comes from the canonical manifold_ node (Embedded: vertex.position, face.{area,
%     normal,centroid}, GlobalVertices/Faces). ASSEMBLED operators (cotan K, mass M, first-order
%     Dirac D) come from operator_ nodes. CONNECTIVITY (face->vertex triples) comes from
%     Surf.Faces. All computations are SEPARATE per hemisphere (the basis is per-hemisphere).
%
% INPUTS:
%     - Data : one of: a results_/matrix filename, a cell of filenames, a matrix [nSrc x nT],
%              or a cell of such matrices. Per-Method input layout (per vertex, both hemispheres):
%                 'gradient'   : scalar field        F [nV  x nT] -> per-FACE tangent vector [3nF x nT] (manifold DEC: #d0)
%                 'divergence' : tangent face field  F [3nF x nT] -> per-VERTEX scalar       [nV  x nT] (manifold DEC: -*0^-1 d0' #' W_F)
%                 'curl'       : tangent face field  F [3nF x nT] -> per-VERTEX vorticity    [nV  x nT] (manifold DEC: -div(N x V))
%                 'laplacian'  : scalar field        F [nV  x nT] -> Laplace-Beltrami        [nV  x nT]
%                 'poisson'    : scalar RHS          F [nV  x nT] -> potential phi            [nV  x nT]  (K phi = M f)
%                 'helmholtz'  : 3-vector field      F [3nV x nT] -> process_helmholtz Hodge decomp (vorticity Curl primary)
%     - OPTIONS : struct (call bst_operators() with no args for defaults), fields below.
%
% OUTPUTS:
%     - OutputFiles : cell-array of created files (or a struct of contents when iTargetStudy='NoSave')
%     - Messages    : error/warning string
%     - isError     : 1 if an error happened
%
% SEE ALSO: bst_eigen, bst_helmholtz, bst_gradient, tess_operators, tess_manifold
%
%     TANGENT STRATUM — two Laplacians, not one. The differential ops (div/curl) use the
%     de Rham / DEC route (delta d + d delta, the Hodge Laplacian). The SMOOTHING / eigen
%     Laplacian for tangent fields is the CONNECTION (Bochner) Laplacian (grad* grad),
%     supplied by the 'Connection Laplacian' operator node and bst_eigen. They differ by
%     Gauss curvature K (Weitzenbock: Delta_Hodge = Delta_connection + Ric, Ric = K g on a
%     surface). Use the de Rham route for div/curl/Helmholtz; the connection route for
%     vector smoothing/spectral analysis.

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

% ===== DEFAULT OPTIONS =====
Def_OPTIONS.Comment      = '';
Def_OPTIONS.Method       = 'gradient';        % {'gradient','divergence','curl','laplacian','poisson','helmholtz'}
Def_OPTIONS.SurfaceFile  = [];                % cortex tying the data to its operators ([] => from the results file)
Def_OPTIONS.Gauge        = 'trivial';         % manifold_ gauge for node resolution
Def_OPTIONS.TimeWindow   = [];                % restrict the input time series
Def_OPTIONS.iTargetStudy = [];                % output study ('NoSave' => return contents, do not save)
Def_OPTIONS.FieldType    = 'auto';            % {'auto','scalar','tangent','ambient'} field stratum

% Return the default options
if (nargin == 0)
    OutputFiles = Def_OPTIONS;
    return;
end
OPTIONS = struct_copy_fields(OPTIONS, Def_OPTIONS, 0);

% ===== PARSE INPUTS =====
OutputFiles = {};
Messages    = [];
isError     = 0;
if ~iscell(Data)
    Data = {Data};
end

bst_progress('start', 'Differential analysis', 'Applying differential operators...', 0, numel(Data));
try
for iData = 1:numel(Data)
    % ----- read the source-mapped data F [nSrc x nT] + its SurfaceFile (mirrors bst_eigen) -----
    isFile = ischar(Data{iData});
    iStudy = [];
    InitFile = [];
    SurfaceFile = OPTIONS.SurfaceFile;
    TimeVector  = [];
    if isFile
        InitFile = Data{iData};
        [~, iStudy, ~, DataType] = bst_get('AnyFile', InitFile);
        if strcmpi(DataType, 'link'), DataType = 'results'; end
        switch DataType
            case 'results'
                R = in_bst_results(InitFile, 1, 'ImageGridAmp', 'Time', 'SurfaceFile');
                F = R.ImageGridAmp;  TimeVector = R.Time;
                if isempty(SurfaceFile), SurfaceFile = R.SurfaceFile; end
            case 'matrix'
                M = in_bst_matrix(InitFile);  F = M.Value;  TimeVector = M.Time;
            otherwise
                Messages = ['bst_operators: unsupported data type: ' DataType];  isError = 1;  break;
        end
        if ~isempty(OPTIONS.TimeWindow) && ~isempty(TimeVector)
            it = bst_closest(OPTIONS.TimeWindow, TimeVector);  it = it(1):it(2);
            TimeVector = TimeVector(it);  F = F(:, it);
        end
    else
        F = Data{iData};
    end
    if isempty(SurfaceFile)
        Messages = 'bst_operators: no SurfaceFile (pass OPTIONS.SurfaceFile for matrix input).';  isError = 1;  break;
    end
    if isempty(F)
        Messages = 'bst_operators: empty source-mapped data; nothing to analyze.';  isError = 1;  break;
    end

    % ----- resolve geometry/operators (per-Method) + dispatch -----
    Surf = in_tess_bst(SurfaceFile, 0);
    nVtot = size(Surf.Vertices, 1);
    nFtot = size(Surf.Faces, 1);
    stratum = i_field_stratum(F, nVtot, nFtot, OPTIONS.FieldType);
    % valid (Method x stratum) pairs; everything else is guarded
    valid = struct('gradient',{{'scalar'}}, 'divergence',{{'tangent','ambient'}}, ...
                   'curl',{{'tangent','ambient'}}, 'laplacian',{{'scalar'}}, ...
                   'poisson',{{'scalar'}}, 'helmholtz',{{'ambient'}});
    m = lower(OPTIONS.Method);
    if isfield(valid, m) && ~ismember(stratum, valid.(m))
        error('bst_operators:badFieldType', ...
            '%s is undefined for a %s field.', OPTIONS.Method, stratum);
    end
    switch lower(OPTIONS.Method)
        case 'gradient'
            Mani  = tess_manifold(SurfaceFile, 'Gauge', OPTIONS.Gauge);   % carries the DEC operators
            f     = i_as_vertex_scalar(F, nVtot, 'gradient');
            [Field, Mag] = bst_gradient(f, Mani);                    % [3nF x nT] per-face native + per-face magnitude
            Result = struct('Method','gradient', 'Field',Field, 'nComponents',3, 'Domain','face', 'Magnitude',Mag);
        case 'laplacian'
            LBO   = tess_operators(SurfaceFile, 'Laplace-Beltrami');
            f     = i_as_vertex_scalar(F, nVtot, 'laplacian');
            Field = i_laplacian(LBO, f, nVtot);                     % [nV x nT]
            Result = struct('Method','laplacian', 'Field',Field, 'nComponents',1);
        case 'poisson'
            LBO   = tess_operators(SurfaceFile, 'Laplace-Beltrami');
            f     = i_as_vertex_scalar(F, nVtot, 'poisson');
            Field = bst_poisson(LBO, f);                            % [nV x nT]
            Result = struct('Method','poisson', 'Field',Field, 'nComponents',1);
        case 'divergence'
            Mani = tess_manifold(SurfaceFile, 'Gauge', OPTIONS.Gauge);
            if strcmp(stratum, 'ambient')
                Dir = tess_operators(SurfaceFile,'Covariant');
                Field = bst_divergence(F, Mani, 'Ambient', Surf, Dir);
            else
                Field = bst_divergence(F, Mani);                   % tangent [3nF]
            end
            Result = struct('Method','divergence', 'Field',Field, 'nComponents',1);
        case 'curl'
            Mani = tess_manifold(SurfaceFile, 'Gauge', OPTIONS.Gauge);
            if strcmp(stratum, 'ambient')
                Dir = tess_operators(SurfaceFile,'Covariant');
                Field = bst_curl(F, Mani, 'Ambient', Surf, Dir);
            else
                Field = bst_curl(F, Mani);                          % tangent [3nF]
            end
            Result = struct('Method','curl', 'Field',Field, 'nComponents',1);
        case 'helmholtz'
            if size(F,1) ~= 3*nVtot
                Messages = sprintf('bst_operators: helmholtz needs a [3nV x nT] vector field (got %d rows, 3nV=%d).', size(F,1), 3*nVtot);
                isError = 1; break;
            end
            Cov = tess_operators(SurfaceFile, 'Covariant');
            H = process_helmholtz('Compute', F, Cov);
            Result = struct('Method','helmholtz', 'Field',H.Curl, 'nComponents',1, 'Helmholtz',H);  % vorticity primary
        otherwise
            Messages = ['bst_operators: unknown Method: ' OPTIONS.Method];  isError = 1;  break;
    end

    % ----- output: save a results_ file, or return contents on NoSave -----
    if isequal(OPTIONS.iTargetStudy, 'NoSave')
        OutputFiles{end+1} = Result; %#ok<AGROW>
    else
        iTarget = OPTIONS.iTargetStudy;
        if isempty(iTarget), iTarget = iStudy; end
        OutFile = i_save_results(iTarget, InitFile, Result, OPTIONS, SurfaceFile, TimeVector);
        if ~isempty(OutFile), OutputFiles{end+1} = OutFile; end %#ok<AGROW>
    end
    bst_progress('inc', 1);
end
catch ME
    bst_progress('stop');
    rethrow(ME);
end
bst_progress('stop');
end


%% ===== infer the field stratum from layout (or honor an explicit FieldType) =====
function stratum = i_field_stratum(F, nVtot, nFtot, explicit)
    if ~strcmpi(explicit, 'auto'), stratum = lower(explicit); return; end
    nr = size(F,1);
    if nr == 3*nVtot
        stratum = 'ambient';
    elseif nr == 3*nFtot
        stratum = 'tangent';
    elseif nr == nVtot
        if nVtot == nFtot
            error('bst_operators:badFieldType', ...
                'Ambiguous layout (nV==nF); set OPTIONS.FieldType explicitly.');
        end
        if ~isreal(F)
            stratum = 'tangent';        % connection (complex) representation
        else
            stratum = 'scalar';
        end
    else
        error('bst_operators:badFieldType', ...
            'Cannot infer field stratum from %d rows (nV=%d, nF=%d); set OPTIONS.FieldType.', nr, nVtot, nFtot);
    end
end


%% ===== coerce input to a per-vertex scalar field [nV x nT] =====
function f = i_as_vertex_scalar(F, nVtot, method)
    if size(F,1) == nVtot
        f = F;
    elseif size(F,1) == 3*nVtot
        % unconstrained source map -> use the vector magnitude as the scalar
        f = sqrt(F(1:3:end,:).^2 + F(2:3:end,:).^2 + F(3:3:end,:).^2);
    else
        error('bst_operators:badInput', ...
            'Method ''%s'' expects a per-vertex scalar [nV=%d x nT] (or [3nV] vector); got %d rows.', ...
            method, nVtot, size(F,1));
    end
end


%% ===== LAPLACIAN: Laplace-Beltrami of a per-vertex scalar, Lf = M^{-1} K f [nV x nT] =====
function Lf = i_laplacian(LBO, f, nVtot)
    Lf = zeros(nVtot, size(f,2));
    for hh = 1:numel(LBO.Operator)
        vH = double(LBO.GlobalVertices{hh}(:));
        K  = LBO.Operator{hh};  M = LBO.Mass{hh};
        Lf(vH,:) = M \ (K * f(vH,:));    % mass-normalized stiffness (cotan K is +semidefinite => this is -Laplace-Beltrami)
    end
end



%% ===== save a results_ source map (mirrors bst_eigen's filter output) =====
function OutFile = i_save_results(iStudy, DataFile, Result, OPTIONS, SurfaceFile, TimeVector)
    OutFile = [];
    if isempty(iStudy), return; end
    if isempty(TimeVector), TimeVector = [0 1]; end
    if size(Result.Field,2) == 1, Result.Field = repmat(Result.Field, 1, 2); end  % viewer needs >=2 time pts
    Comment = OPTIONS.Comment;
    if isempty(Comment), Comment = ['Differential | ' Result.Method]; end
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImageGridAmp  = Result.Field;
    ResultsMat.Time          = TimeVector(:)';
    ResultsMat.nComponents   = Result.nComponents;
    ResultsMat.Comment       = Comment;
    ResultsMat.SurfaceFile   = SurfaceFile;
    ResultsMat.Function      = 'bst_operators';
    ResultsMat.DataFile      = DataFile;
    ResultsMat.Options       = struct('Method', Result.Method);
    if isfield(Result, 'Domain'), ResultsMat.Options.Domain = Result.Domain; end   % 'face' for gradient
    sStudy  = bst_get('Study', iStudy);
    OutFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['results_differential_' Result.Method '_']);
    bst_save(OutFile, ResultsMat, 'v6');
    db_add_data(iStudy, OutFile, ResultsMat);
end
