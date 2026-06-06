function CompHM = bst_face_eigenmode_leadfield(FaceHeadModelFile, varargin)
% BST_FACE_EIGENMODE_LEADFIELD  Compose a face-based leadfield onto scalar LBO eigenmodes.
%
% USAGE:
%   CompHM = bst_face_eigenmode_leadfield(FaceHeadModelFile)
%   CompHM = bst_face_eigenmode_leadfield(FaceHeadModelFile, 'nModes', 200)
%   CompHM = bst_face_eigenmode_leadfield(FaceHeadModelFile, 'Register', true)
%
% DESCRIPTION:
%   Face-based analogue of bst_eigenmode_leadfield.  Takes a face headmodel
%   (from bst_face_headmodel), projects its leadfield onto LBO eigenmodes
%   defined on the face mesh, and saves a composed headmodel compatible with
%   the full Brainstorm eigenmode inverse pipeline.
%
%   The composed leadfield is:
%       L̃[:, k] = L_face * Ψ_k = Σ_f Ψ_k(f) · G(x_f) · n̂_f · A_f
%
%   Each column is the sensor pattern of a distributed current sheet whose
%   density follows the k-th spatial eigenmode — a proper quadrature of
%   the surface integral with exact face normals and area weighting.
%
%   This is geometrically complete: the forward model encodes the exact
%   face normal (from the trivial-connection frame), the face area, and
%   the smoothness of the eigenmode basis.
%
% EIGENMODES:
%   Vertex LBO eigenmodes from the surface (bst_eigenmodes_ensure) are
%   interpolated to faces by corner averaging and M_f-orthonormalized.
%   For reconstruction, they are also interpolated back to vertices so
%   bst_eigenmode_reconstruct can produce a vertex-indexed imaging kernel
%   that works with the standard Brainstorm display pipeline.
%
% INPUTS:
%   FaceHeadModelFile  Brainstorm relative path to a face headmodel
%                      (produced by bst_face_headmodel, isFaceBased=1)
%
% OPTIONS (name-value):
%   'nModes'    K, number of eigenmodes to retain (default: 200)
%   'Register'  logical, register in DB and set as current HM (default: true)
%
% OUTPUT:
%   CompHM  composed head-model struct (also saved if Register=true)
%           Key fields:
%             .Gain         [nCh x K]    composed leadfield L̃ = L_face * Ψ
%             .Eigenvalues  [K x 1]      LBO eigenvalues
%             .nModes       K
%             .ModeIndices  [K x 1]      which eigenmodes were selected
%             .isEigenmode  1            flag for bst_inverse_eigenmodes
%             .isFaceBased  1            flag for bst_eigenmode_reconstruct
%             .PhiVertices  [nV x K]     Ψ interpolated to vertices
%                                        (used by bst_eigenmode_reconstruct
%                                         when isFaceBased=1)
%
% PIPELINE:
%   bst_face_headmodel          → face headmodel  (L_face [nCh x nF])
%   bst_face_eigenmode_leadfield → composed HM    (L̃ = L_face·Ψ [nCh x K])
%   bst_inverse_eigenmodes       → mode kernel    (θ-space [K x nCh])
%   bst_eigenmode_reconstruct    → imaging kernel ([nV x nCh])
%   [apply to data]              → s(v,t)          vertex source maps
%
% SEE ALSO: bst_face_headmodel, bst_face_leadfield, bst_eigenmode_leadfield,
%           bst_inverse_eigenmodes, bst_eigenmode_reconstruct
%           dev/references/face_based_source_model.md
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
nModes   = 200;
Register = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'nmodes',   nModes   = varargin{k+1};
        case 'register', Register = logical(varargin{k+1});
    end
end

%% ── Load face headmodel ──────────────────────────────────────────────────
FHM = in_bst_headmodel(FaceHeadModelFile, 0);
if ~isfield(FHM, 'isFaceBased') || ~FHM.isFaceBased
    error('bst_face_eigenmode_leadfield:NotFaceBased', ...
        'Input must be a face-based headmodel (isFaceBased=1). Use bst_face_headmodel first.');
end

L_face      = double(FHM.Gain);           % [nCh x nF]
SurfaceFile = FHM.SurfaceFile;
nCh         = size(L_face, 1);
nF          = size(L_face, 2);

fprintf('bst_face_eigenmode_leadfield: %d channels, %d faces, surface %s\n', ...
    nCh, nF, SurfaceFile);

%% ── Surface and mesh structures ─────────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile);
Faces   = double(TessMat.Faces);
nV      = size(TessMat.Vertices, 1);

[~, lH_v] = tess_hemisplit(TessMat);  lH_v = lH_v(:);
inL_v      = false(nV,1); inL_v(lH_v) = true;
lH_f       = find(all(inL_v(Faces), 2));   % LH face indices
nLHF       = numel(lH_f);

%% ── Vertex LBO eigenmodes → face eigenmodes ──────────────────────────────
% 1. Ensure vertex LBO modes are computed on this surface
Eig = bst_eigenmodes_ensure(SurfaceFile, nModes);
K   = min(nModes, Eig.nModes);

% Canonical global order (ascending eigenvalue)
if isfield(Eig,'Order') && ~isempty(Eig.Order)
    order = double(Eig.Order(:));
else
    [~, order] = sort(double(Eig.Values(:)), 'ascend');
end
sel        = order(1:K);
Phi_v      = double(Eig.Vectors(:, sel));   % [nV x K]  vertex eigenmodes
lam        = double(Eig.Values(sel));        % [K x 1]   eigenvalues

% 2. Interpolate vertex modes to LH faces by corner average
lh_vmap  = zeros(nV,1); lh_vmap(lH_v) = 1:numel(lH_v);
FacesLH  = Faces(lH_f,:);                   % [nLHF x 3]

% Map: vertex indices in LH-local space
iv1 = lh_vmap(FacesLH(:,1));
iv2 = lh_vmap(FacesLH(:,2));
iv3 = lh_vmap(FacesLH(:,3));

Phi_v_lh   = Phi_v(lH_v, :);               % [nLH x K]
Phi_f_raw  = (Phi_v_lh(iv1,:) + Phi_v_lh(iv2,:) + Phi_v_lh(iv3,:)) / 3;  % [nLHF x K]

% 3. M_f-orthonormalize under face area matrix
% Load DEC operators for face areas
clear mex
h_dec = nxr_compute('create', TessMat.Vertices, Faces);
dec   = nxr_compute('assembleDECOperators', h_dec);
nxr_compute('destroy', h_dec);
M_f_lh = dec.hodge2(lH_f, lH_f);           % [nLHF x nLHF] diagonal face areas

G_gram = Phi_f_raw' * M_f_lh * Phi_f_raw;  % [K x K] Gram matrix
[R, p] = chol(G_gram);
if p ~= 0
    warning('bst_face_eigenmode_leadfield:GramNotPD', ...
        'Face Gram matrix not PD (condition issue); using mass-normalization only.');
    nrm = sqrt(diag(G_gram))';
    Phi_f = Phi_f_raw ./ nrm;
else
    Phi_f = Phi_f_raw / R;                  % [nLHF x K]  M_f-orthonormal
end

fprintf('  Face eigenmodes: %d  (λ range [%.0f … %.0f])\n', K, min(lam), max(lam));

%% ── Composed leadfield L̃ = L_face · Ψ ────────────────────────────────────
L_tilde = L_face(:, lH_f) * Phi_f;         % [nCh x K]
fprintf('  L̃ = L_face · Ψ:  [%d x %d]\n', size(L_tilde,1), size(L_tilde,2));

%% ── Face → vertex interpolation for reconstruction ───────────────────────
% bst_eigenmode_reconstruct expects Phi [nV x K] to produce a vertex-indexed
% imaging kernel.  We interpolate Phi_f back to vertices by averaging over
% all incident LH faces; RH vertices receive zeros.
% This is the "display interpolation" step — the forward model used exact
% face geometry; reconstruction maps back to vertices for the viewer.
Phi_V = zeros(nV, K);
for vi = 1:numel(lH_v)
    v    = lH_v(vi);
    % Find LH faces that contain this vertex (in local lH_f indexing)
    fmask = (FacesLH(:,1) == v) | (FacesLH(:,2) == v) | (FacesLH(:,3) == v);
    if any(fmask)
        Phi_V(v,:) = mean(Phi_f(fmask,:), 1);
    end
end

%% ── Build composed head-model struct ─────────────────────────────────────
CompHM = struct();
CompHM.MEGMethod     = FHM.MEGMethod;
CompHM.EEGMethod     = '';
CompHM.ECOGMethod    = '';
CompHM.SEEGMethod    = '';
CompHM.NIRSMethod    = '';
CompHM.Gain          = L_tilde;
CompHM.Comment       = sprintf('Face Eigenmode leadfield (%d modes) | %s', K, SurfaceFile);
CompHM.HeadModelType = 'surface';
CompHM.GridLoc       = FHM.GridLoc;         % face centroids
CompHM.GridOrient    = FHM.GridOrient;      % face normals
CompHM.GridAtlas     = [];
CompHM.GridOptions   = [];
CompHM.SurfaceFile   = file_win2unix(SurfaceFile);
CompHM.Param         = FHM.Param;

% Eigenmode fields (required by bst_inverse_eigenmodes)
CompHM.isEigenmode   = 1;
CompHM.nModes        = K;
CompHM.Eigenvalues   = lam(:);
CompHM.ModeIndices   = sel(:);

% Face-based fields
CompHM.isFaceBased   = 1;
CompHM.GridAreas     = FHM.GridAreas;       % face areas [m²]
CompHM.GridU         = FHM.GridU;           % trivial-connection e1
CompHM.GridV         = FHM.GridV;           % trivial-connection e2
CompHM.PhiVertices   = Phi_V;              % [nV x K]  for bst_eigenmode_reconstruct

CompHM = bst_history('add', CompHM, 'compute', ...
    sprintf('Face eigenmode leadfield | %s | %d modes', FaceHeadModelFile, K));

%% ── Save and register ────────────────────────────────────────────────────
if ~Register
    return
end

ProtocolInfo   = bst_get('ProtocolInfo');
[sStudy, iStudy] = bst_get('HeadModelFile', FaceHeadModelFile);

OutFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
    'headmodel_face_eigenmode');
if ~strcmpi(OutFile(end-3:end), '.mat'), OutFile = [OutFile '.mat']; end

bst_save(OutFile, CompHM, 'v7');
fprintf('  Saved → %s\n', OutFile);

newHM              = db_template('HeadModel');
newHM.FileName     = file_win2unix(strrep(OutFile, ProtocolInfo.STUDIES, ''));
newHM.Comment      = CompHM.Comment;
newHM.HeadModelType = 'surface';
newHM.MEGMethod    = CompHM.MEGMethod;
newHM.EEGMethod    = ''; newHM.ECOGMethod = ''; newHM.SEEGMethod = ''; newHM.NIRSMethod = '';

iNew = length(sStudy.HeadModel) + 1;
sStudy.HeadModel(iNew) = newHM;
sStudy.iHeadModel      = iNew;
bst_set('Study', iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', iStudy);

fprintf('  Registered as study head model #%d\n', iNew);
CompHM.FileName = newHM.FileName;

end
