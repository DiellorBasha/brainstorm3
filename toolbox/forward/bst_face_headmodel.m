function HeadModelFile = bst_face_headmodel(BaseHeadModelFile, varargin)
% BST_FACE_HEADMODEL  Compute and register a face-based MEG head model.
%
% USAGE:
%   HeadModelFile = bst_face_headmodel(BaseHeadModelFile)
%   HeadModelFile = bst_face_headmodel(BaseHeadModelFile, 'Mode', 'constrained')
%   HeadModelFile = bst_face_headmodel(BaseHeadModelFile, 'Mode', 'loose')
%   HeadModelFile = bst_face_headmodel(BaseHeadModelFile, 'BlockSize', 500)
%
% DESCRIPTION:
%   Computes a face-based constrained (or loose) MEG leadfield and registers
%   it as a new head model in the Brainstorm database.
%
%   The source model uses triangular faces as source elements rather than
%   vertices.  Each face contributes one (constrained) or three (loose)
%   leadfield columns.  The tangent frame {n̂_f, U_f, V_f} comes entirely
%   from tess_tangents (FreeSurfer trivial-connection), ensuring:
%     • n̂_f is the exact face normal — no vertex averaging
%     • U_f, V_f are smooth, globally consistent — correct gauge for
%       wave analysis, DEC gradient, and connection-Laplacian eigenmodes
%
%   The sphere parameters (Param) and channel structure are reused from the
%   base os_meg headmodel — no re-fitting of spheres required.
%
% INPUTS:
%   BaseHeadModelFile  Brainstorm relative path to an existing os_meg surface
%                      headmodel for the same subject/condition.  Used to
%                      obtain sphere parameters, channel mapping, and the
%                      cortical surface file.
%
% OPTIONS (name-value):
%   'Mode'       'constrained' (default) | 'loose' | 'unconstrained'
%                unconstrained: 3 Cartesian columns/face = raw Sarvas at centroid
%                               [nCh x 3F], full-3D, geometry from tess_manifold
%                               (the rigorous full-unconstrained inverse path)
%                constrained:   1 column/face along n̂_f (normal current flux)
%                loose:         3 columns/face [n̂_f*A_f, U_f, V_f] (legacy frame)
%   'BlockSize'  faces per Sarvas batch call (default 500)
%
% OUTPUT:
%   HeadModelFile  Brainstorm relative path of the saved file, registered
%                  in the database and set as the current head model.
%
% SAVED FILE FIELDS:
%   Standard Brainstorm headmodel fields:
%     Gain, GridLoc, GridOrient, HeadModelType='surface',
%     MEGMethod, Comment, SurfaceFile, Param, History
%   Additional face-based fields:
%     GridAreas   [nF x 1]  face areas [m²]
%     GridU       [nF x 3]  trivial-connection tangent e1 (smooth gauge)
%     GridV       [nF x 3]  trivial-connection tangent e2 = n̂_f × U_f
%     isFaceBased  1         flag for downstream face-aware routines
%     nComponents  1 (constrained) or 3 (loose)
%
% SEE ALSO: bst_face_leadfield, tess_tangents, bst_eigenmode_leadfield
%           dev/references/face_based_source_model.md
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
Mode      = 'constrained';
BlockSize = 500;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'mode',      Mode      = lower(varargin{k+1});
        case 'blocksize', BlockSize = varargin{k+1};
    end
end

%% ── Load base head model ─────────────────────────────────────────────────
% The base os_meg headmodel provides the sphere parameters (Param) and
% channel mapping that we reuse — no re-fitting required.
BaseHM = in_bst_headmodel(BaseHeadModelFile, 0);

if ~strcmpi(BaseHM.MEGMethod, 'os_meg') && ~strcmpi(BaseHM.MEGMethod, 'meg_sphere')
    error('bst_face_headmodel:MEGMethod', ...
        'Base headmodel must use os_meg or meg_sphere (got ''%s'').', BaseHM.MEGMethod);
end
if ~strcmpi(BaseHM.HeadModelType, 'surface')
    error('bst_face_headmodel:HeadModelType', ...
        'Face-based model requires a surface-type base headmodel.');
end

SurfaceFile = BaseHM.SurfaceFile;

%% ── Channel structure and sphere parameters ──────────────────────────────
% Get study to find channel file
[sStudy, iStudy] = bst_get('HeadModelFile', BaseHeadModelFile);
if isempty(sStudy)
    error('bst_face_headmodel:StudyNotFound', ...
        'Could not locate study for headmodel: %s', BaseHeadModelFile);
end

ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
ChannelMat  = in_bst_channel(ChannelFile);

% Select MEG channels — type 'MEG' (gradiometers/magnetometers, not REF).
% The base Gain may have more rows (all 340 channels including EEG/EOG/REF
% as NaN rows); we only need the MEG subset for the face forward model.
% BaseHM.Param has one entry per channel in the full channel file, so we
% index it with iMEG directly.
iMEG = find(strcmpi({ChannelMat.Channel.Type}, 'MEG'));
if isempty(iMEG)
    error('bst_face_headmodel:NoMEG', 'No MEG channels found in channel file.');
end
if numel(BaseHM.Param) < max(iMEG)
    error('bst_face_headmodel:ParamMismatch', ...
        'BaseHM.Param has %d entries but channel index %d requested.', ...
        numel(BaseHM.Param), max(iMEG));
end

Channel_meg = ChannelMat.Channel(iMEG);
Param_meg   = BaseHM.Param(iMEG);

fprintf('bst_face_headmodel: %d MEG channels, surface %s, mode=%s\n', ...
    numel(Channel_meg), SurfaceFile, Mode);

%% ── Compute face-based leadfield ─────────────────────────────────────────

[L_face, FaceGeom] = bst_face_leadfield(SurfaceFile, Channel_meg, Param_meg, ...
    'Mode', Mode, 'BlockSize', BlockSize);

nF    = size(FaceGeom.Centroids, 1);
nComp = 1 + 2*(strcmpi(Mode, 'loose') || strcmpi(Mode, 'unconstrained'));   % 1 constrained, 3 loose/unconstrained

fprintf('bst_face_headmodel: Gain = [%d x %d]  (%d faces, %d component(s))\n', ...
    size(L_face,1), size(L_face,2), nF, nComp);

%% ── Build and save headmodel struct ─────────────────────────────────────

ProtocolInfo = bst_get('ProtocolInfo');

% Output filename
if     strcmpi(Mode,'loose'),         modeTag = 'loose';
elseif strcmpi(Mode,'unconstrained'), modeTag = 'unconstrained';
else,                                 modeTag = 'constrained';
end
OutFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
    sprintf('headmodel_face_%s_%s', BaseHM.MEGMethod, modeTag));
% GetNewFilename already appends .mat — do not add again
if ~strcmpi(OutFile(end-3:end), '.mat')
    OutFile = [OutFile, '.mat'];
end

% Standard Brainstorm headmodel fields
HM = struct();
HM.MEGMethod     = BaseHM.MEGMethod;
HM.EEGMethod     = '';
HM.ECOGMethod    = '';
HM.SEEGMethod    = '';
HM.NIRSMethod    = '';
HM.Gain          = L_face;
HM.Comment       = sprintf('Face-%s (%s) | %d faces | %s', ...
    modeTag, BaseHM.MEGMethod, nF, SurfaceFile);
HM.HeadModelType = 'surface';
HM.GridLoc       = FaceGeom.Centroids;
HM.GridOrient    = FaceGeom.Normals;
HM.GridAtlas     = [];
HM.GridOptions   = [];
HM.SurfaceFile   = file_win2unix(SurfaceFile);
HM.Param         = Param_meg;

% Additional face-based fields (not present in vertex headmodels)
HM.GridAreas     = FaceGeom.Areas;    % [nF x 1]  face areas [m²]
HM.GridU         = []; if isfield(FaceGeom,'U'), HM.GridU = FaceGeom.U; end   % [nF x 3] frame e1 (legacy modes only)
HM.GridV         = []; if isfield(FaceGeom,'V'), HM.GridV = FaceGeom.V; end   % [nF x 3] frame e2 (legacy modes only)
HM.isFaceBased   = 1;                 % flag: face-indexed sources, not vertex-indexed
HM.nComponents   = nComp;

HM = bst_history('add', HM, 'compute', ...
    sprintf('Face-based %s headmodel from %s | %s mode | %d faces', ...
        BaseHM.MEGMethod, BaseHeadModelFile, modeTag, nF));

bst_save(OutFile, HM, 'v7');
fprintf('bst_face_headmodel: saved → %s\n', OutFile);

%% ── Register in Brainstorm database ─────────────────────────────────────

newHM              = db_template('HeadModel');
newHM.FileName     = file_win2unix(strrep(OutFile, ProtocolInfo.STUDIES, ''));
newHM.Comment      = HM.Comment;
newHM.HeadModelType = 'surface';
newHM.MEGMethod    = HM.MEGMethod;
newHM.EEGMethod    = '';
newHM.ECOGMethod   = '';
newHM.SEEGMethod   = '';
newHM.NIRSMethod   = '';

iNew = length(sStudy.HeadModel) + 1;
sStudy.HeadModel(iNew) = newHM;
sStudy.iHeadModel      = iNew;

bst_set('Study', iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', iStudy);

HeadModelFile = newHM.FileName;
fprintf('bst_face_headmodel: registered as study head model #%d\n', iNew);

end


function v = iff(cond, a, b)
if cond, v = a; else, v = b; end
end
