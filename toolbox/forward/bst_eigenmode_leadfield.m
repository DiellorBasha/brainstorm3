function CompHM = bst_eigenmode_leadfield(HeadModel, Eigenmodes, varargin)
% BST_EIGENMODE_LEADFIELD: Compose a leadfield into the LBO eigenmode basis.
%
% USAGE:  CompHM = bst_eigenmode_leadfield(HeadModel, Eigenmodes, 'nModes', K)
%
% DESCRIPTION:
%     Forward solution for eigenmode-space source mapping (GBF). Takes a base
%     surface head model and the precomputed surface eigenmodes and returns a
%     composed head-model struct whose Gain is the eigenmode leadfield
%         L~ = L * Phi            [nChannels x K]
%     where L is the constrained (surface-normal) leadfield and Phi = Eigenmodes
%     truncated to K modes. Each column of L~ is the sensor topography of one
%     eigenmode. This is strictly a forward operation; the inverse is separate.
%
% INPUTS:
%     HeadModel  : base head-model struct (from in_bst_headmodel, ApplyOrient=0):
%                  .Gain [nCh x 3*nVert] unconstrained, .GridOrient [nVert x 3],
%                  .GridAtlas, .HeadModelType, .SurfaceFile, .Comment
%     Eigenmodes : struct from in_tess_eigenmodes: .Vectors [nVert x nModes],
%                  .Values [nModes x 1], .nModes
% OPTIONS:
%     'nModes' : number of leading modes to keep (default: all; clamped to available)
%
% OUTPUT:
%     CompHM : composed head-model struct ready to save as headmodel_eigenmode_*.mat
%
% Authors: Diellor Basha, 2026

nModes = [];
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'nmodes', nModes = varargin{i+1};
    end
end

% Constrained leadfield: [nCh x nVert]
Lc = bst_gain_orient(double(HeadModel.Gain), HeadModel.GridOrient, ...
    getfield_default(HeadModel, 'GridAtlas', []));

Phi    = double(Eigenmodes.Vectors);     % [nVert x nModesAll]
Values = double(Eigenmodes.Values(:));   % [nModesAll x 1]
nVert  = size(Phi, 1);
nModesAll = size(Phi, 2);

if size(Lc,2) ~= nVert
    error('bst_eigenmode_leadfield:VertexMismatch', ...
        'Leadfield has %d sources but eigenmodes have %d vertices.', size(Lc,2), nVert);
end

% Clamp K
if isempty(nModes) || nModes <= 0
    K = nModesAll;
else
    K = min(nModes, nModesAll);
end
Phi = Phi(:, 1:K);

% Compose: L~ = L * Phi   [nCh x K]
L_tilde = Lc * Phi;

% Build composed head-model struct
CompHM = HeadModel;
CompHM.Gain        = L_tilde;
CompHM.GridLoc     = [];
CompHM.GridOrient  = [];
CompHM.GridAtlas   = [];
CompHM.isEigenmode = 1;
CompHM.nModes      = K;
CompHM.Eigenvalues = Values(1:K);
CompHM.HeadModelType = 'surface';
CompHM.Comment     = sprintf('Eigenmode leadfield (%d modes) | %s', K, ...
    getfield_default(HeadModel, 'Comment', ''));
end

function v = getfield_default(s, f, d)
if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
