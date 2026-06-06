function [Z_joint, phi_F, phi_T, confidence] = bst_fiedler_joint_field(Theta, FaceIndices, SurfaceFile, s_alpha, varargin)
% BST_FIEDLER_JOINT_FIELD  Form the joint Fiedler torus field from eigenmode CWT.
%
% USAGE:
%   [Z, phi_F, phi_T, conf] = bst_fiedler_joint_field(Theta, FaceIndices, SurfaceFile, s_alpha)
%
% DESCRIPTION:
%   Constructs the joint Fiedler field Z_joint(x,t) = conj(u₁(x)) · θ₁(s,t)
%   where:
%     u₁(x) ∈ ℂ  — face Fiedler eigenmode from TessMat.nxr (spatial Fiedler)
%     θ₁(s,t) ∈ ℂ — CWT coefficient of the DOMINANT FACE EIGENMODE at scale s
%
%   arg(Z_joint(x,t)) = arg(θ₁(s,t)) − arg(u₁(x))
%                     ≈ ω₀t − k·φ_F(x) + const   [wave phase in Fiedler coordinates]
%
%   where φ_F(x) = arg(u₁(x)) is the Fiedler longitude — a smooth scalar
%   coordinate on the cortex that is continuous across sulcal walls.
%
%   Sign ambiguity is automatic: conj(u₁(x)) demodulates the wave relative
%   to the Fiedler frame, absorbing sulcal π-jumps in u₁ not in the wave.
%
% INPUTS:
%   Theta       [K × nScales × nTime] complex — from bst_eigenmode_cwt_inverse
%               OR [K × nTime] complex — eigenmode analytic signal (Hilbert)
%   FaceIndices [nLHF × 1]   global face indices from HeadModel.FaceIndices
%   SurfaceFile Brainstorm cortex surface file path
%   s_alpha     integer — CWT scale index to use (from bst_dispersion_fit)
%               ignored if Theta is [K × nTime]
%
% OPTIONS (name-value):
%   'ModeIndex'  which eigenmode to use as temporal Fiedler (default: 1 = dominant)
%   'Verbose'    logical (default: true)
%
% OUTPUTS:
%   Z_joint    [nLHF × nTime] complex — joint Fiedler field
%   phi_F      [nLHF × 1]  real — Fiedler longitude = arg(u₁(x)) ∈ (−π,π]
%   phi_T      [1 × nTime] real — temporal phase = arg(θ₁(s,t))
%   confidence [nLHF × nTime] real — |u₁(x)| · |θ₁(s,t)| product
%
% AUTHORS: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
ModeIndex = 1;
Verbose   = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'modeindex', ModeIndex = varargin{k+1};
        case 'verbose',   Verbose   = logical(varargin{k+1});
    end
end

%% ── Load face Fiedler eigenvector from TessMat.nxr ───────────────────────
TessMat = in_tess_bst(SurfaceFile, 0);
if ~isfield(TessMat,'nxr') || isempty(TessMat.nxr)
    error('bst_fiedler_joint_field:NoNxr', ...
        'TessMat.nxr missing. Run tess_nxr_populate first.');
end
NxrData = TessMat.nxr;

% Map FaceIndices to LH-local index in ConnEigLH
lH_f_nxr = NxrData.lH_f;
[~, lhf_loc] = ismember(FaceIndices(:), lH_f_nxr);
if any(lhf_loc == 0)
    warning('bst_fiedler_joint_field: %d FaceIndices not in nxr LH face set.', sum(lhf_loc==0));
end

nLHF = numel(FaceIndices);
u1_full = zeros(nLHF, 1, 'like', 1i);
valid_f = lhf_loc > 0;
u1_full(valid_f) = NxrData.ConnEigLH.Vectors(lhf_loc(valid_f), ModeIndex);

% Fiedler longitude: smooth scalar coordinate on cortex
phi_F = angle(u1_full);      % [nLHF × 1] ∈ (−π, π]
amp_F = abs(u1_full);        % confidence per face

%% ── Extract temporal Fiedler: θ₁(s,t) at selected scale ─────────────────
ndims_Theta = ndims(Theta);
if ndims_Theta == 3
    % Theta [K × nScales × nTime]
    theta1 = squeeze(Theta(ModeIndex, s_alpha, :)).';   % [1 × nTime]
else
    % Theta [K × nTime] — from Hilbert inverse (no scale dimension)
    theta1 = Theta(ModeIndex, :);                        % [1 × nTime]
end

phi_T  = angle(theta1);   % [1 × nTime] temporal phase
amp_T  = abs(theta1);     % [1 × nTime] burst amplitude

%% ── Joint field: outer product in complex notation ───────────────────────
% Z_joint(x,t) = conj(u₁(x)) · θ₁(s,t)
% arg(Z_joint) = arg(θ₁) − arg(u₁) ≈ wave phase in Fiedler coordinates
Z_joint    = conj(u1_full) .* theta1;          % [nLHF × nTime] broadcast
confidence = amp_F .* amp_T;                   % [nLHF × nTime] joint confidence

if Verbose
    fprintf('bst_fiedler_joint_field: %d LH faces  mode=%d\n', nLHF, ModeIndex);
    fprintf('  Fiedler |u₁| range: [%.3f %.3f]  phi_F range: [%.2f %.2f] rad\n', ...
        min(amp_F), max(amp_F), min(phi_F), max(phi_F));
    fprintf('  |θ₁| mean=%.3e  phi_T range: [%.2f %.2f] rad\n', ...
        mean(amp_T), min(phi_T), max(phi_T));
end
end
