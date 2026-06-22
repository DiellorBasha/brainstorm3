function [Phi, lam, info] = bst_helmet_eigenmodes(ChannelFile, varargin)
% BST_HELMET_EIGENMODES  LBO eigenmodes of the MEG sensor helmet surface.
%
% USAGE:
%   [Phi, lam, info] = bst_helmet_eigenmodes(ChannelFile)
%   [Phi, lam, info] = bst_helmet_eigenmodes(ChannelFile, 'nModes', 60)
%
% DESCRIPTION:
%   Computes the Laplace-Beltrami eigenmodes of the MEG sensor array
%   surface — the continuous analogue of spherical harmonics adapted to the
%   actual helmet geometry.  The sensor positions define a point cloud on a
%   roughly hemi-ellipsoidal surface; a convex-hull triangulation is used as
%   the manifold approximation.
%
%   The resulting eigenmodes are smooth spatial patterns on the sensor array,
%   ordered by spatial frequency (λ_1 ≤ λ_2 ≤ …).  Low-j modes span large-
%   scale, spatially coherent activations (brain signal); high-j modes
%   resolve sensor-scale fluctuations (incoherent noise).
%
%   This is the discrete, geometry-adapted generalisation of the Signal
%   Space Separation (SSS) internal basis: for a perfectly spherical helmet
%   the eigenmodes converge to spherical harmonics and the projection onto
%   the first J modes is equivalent to SSS.
%
% INPUT:
%   ChannelFile  Brainstorm relative path to a channel file
%
% OPTIONS (name-value):
%   'nModes'    J, number of eigenmodes (default: 60)
%   'ChanType'  sensor type string (default: 'MEG')
%   'Verbose'   logical (default: true)
%
% OUTPUTS:
%   Phi   [nCh x J]  eigenvectors at sensor positions (M-normalised)
%   lam   [J x 1]    eigenvalues (ascending)
%   info  struct:
%           .nCh        — number of MEG channels used
%           .nModes     — J
%           .SensorPos  — [nCh x 3] sensor positions [m]
%           .Faces      — convex-hull triangulation
%           .lam_ratio  — lam(J)/lam(1) — dynamic range of spatial frequencies
%
% SPATIAL FILTERING:
%   Project sensor data onto the first J smooth modes:
%       d_filt = Phi * (Phi' * d)       [nCh x nTime]
%   This removes spatially incoherent noise (high-j modes) while preserving
%   spatially smooth brain signals.  J=60 retains ~22% of sensor modes;
%   typical SSS uses ~80 of 306 (26%).
%
% SEE ALSO: bst_eigenmode_analytic_inverse, nxr_compute
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
nModes   = 60;
ChanType = 'MEG';
Verbose  = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'nmodes',   nModes   = varargin{k+1};
        case 'chantype', ChanType = varargin{k+1};
        case 'verbose',  Verbose  = varargin{k+1};
    end
end

%% ── Load channel file and extract sensor positions ───────────────────────
ChanMat = in_bst_channel(ChannelFile);
iCh     = find(strcmpi({ChanMat.Channel.Type}, ChanType));
if isempty(iCh)
    error('bst_helmet_eigenmodes:NoChan', ...
        'No channels of type ''%s'' found in %s.', ChanType, ChannelFile);
end
nCh = numel(iCh);

% Use first coil position for each sensor (column 1 of .Loc [3 x nCoils])
pos = zeros(nCh, 3);
for k = 1:nCh
    loc = ChanMat.Channel(iCh(k)).Loc;
    pos(k, :) = loc(1:3, 1)';
end

if Verbose
    fprintf('bst_helmet_eigenmodes: %d %s channels  (pos range [%.0f %.0f] mm)\n', ...
        nCh, ChanType, min(vecnorm(pos,2,2))*1000, max(vecnorm(pos,2,2))*1000);
end

%% ── Build helmet surface mesh via convex hull ────────────────────────────
% Convex hull may leave some sensors unreferenced (not on hull surface).
% Strategy: compute hull, reindex to hull-only vertices, compute LBO there,
% then interpolate back to all nCh sensors via nearest-hull-vertex mapping.
[F_raw, ~] = convhull(pos(:,1), pos(:,2), pos(:,3));
hull_v  = unique(F_raw(:));          % indices into pos (1-based)
nHull   = numel(hull_v);
pos_hull = pos(hull_v, :);           % [nHull x 3]

% Reindex faces to hull-local indices
fmap = zeros(nCh, 1);
fmap(hull_v) = 1:nHull;
F_hull = fmap(F_raw);                % [nFaces x 3] hull-local

if Verbose
    fprintf('  Convex hull: %d faces, %d/%d sensors on hull\n', ...
        size(F_hull,1), nHull, nCh);
end

%% ── Compute LBO on hull vertices via nxr-compute ─────────────────────────
bst_plugin('Load', 'nxr-compute');
clear mex
h_dec = nxr_compute('create', pos_hull, F_hull);
dec   = nxr_compute('assembleDECOperators', h_dec);
nxr_compute('destroy', h_dec);

L_helm = dec.d0' * dec.hodge1 * dec.d0;   % stiffness [nHull x nHull]
M_helm = dec.hodge0;                        % mass      [nHull x nHull]

%% ── Eigendecomposition on hull vertices ──────────────────────────────────
nModes = min(nModes, nHull - 1);
[V, D] = eigs(L_helm, M_helm, nModes + 1, 'smallestabs');
lam_all = real(diag(D));
[lam_all, ord] = sort(lam_all, 'ascend');
V = real(V(:, ord));

% Skip DC mode (λ ≈ 0)
dc = find(lam_all < 1e-6 * max(abs(lam_all)), 1, 'last');
if isempty(dc), dc = 1; end
sel   = (dc+1) : min(dc+nModes, numel(lam_all));
Phi_hull = V(:, sel);    % [nHull x J]
lam      = lam_all(sel); % [J x 1]

%% ── Interpolate to all nCh sensors via nearest-hull-vertex ───────────────
% For sensors on the hull: direct assignment.
% For off-hull sensors: nearest hull vertex in 3-D.
Phi = zeros(nCh, numel(lam));
Phi(hull_v, :) = Phi_hull;

off_hull = setdiff(1:nCh, hull_v);
if ~isempty(off_hull)
    % Nearest hull-vertex interpolation
    dists = pdist2(pos(off_hull,:), pos_hull);   % [nOff x nHull]
    [~, nn] = min(dists, [], 2);
    Phi(off_hull, :) = Phi_hull(nn, :);
    if Verbose
        fprintf('  %d off-hull sensors interpolated via nearest hull vertex\n', ...
            numel(off_hull));
    end
end

if Verbose
    fprintf('  Helmet LBO: %d modes  λ∈[%.1f … %.1f]  ratio=%.0f×\n', ...
        numel(lam), min(lam), max(lam), max(lam)/max(min(lam),eps));
end

%% ── Pack info ────────────────────────────────────────────────────────────
info.nCh       = nCh;
info.nModes    = numel(lam);
info.SensorPos = pos;
info.Faces     = F_hull;
info.lam_ratio = max(lam) / max(min(lam), eps);
end
