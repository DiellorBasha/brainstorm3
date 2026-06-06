function [Z_joint, phi_F, phi_T, confidence] = bst_fiedler_joint_field(s_face, FaceIndices, SurfaceFile, varargin)
% BST_FIEDLER_JOINT_FIELD  Demodulate a face source field into the Fiedler gauge.
%
% USAGE:
%   [Z, phi_F, phi_T, conf] = bst_fiedler_joint_field(s_face, FaceIndices, SurfaceFile)
%
% DESCRIPTION:
%   Forms the joint Fiedler field by demodulating the FULL reconstructed face
%   source field against the unit spatial-Fiedler phasor:
%
%       Z_joint(x,t) = conj( u₁(x)/|u₁(x)| ) · s_face(x,t)
%
%   where:
%     s_face(x,t) — full face 2-form reconstruction (all modes superposed),
%                   e.g. FaceAmp(:, s_alpha, :) from bst_eigenmode_cwt_inverse
%                   or FaceGridAmp from bst_eigenmode_analytic_inverse
%     u₁(x) ∈ ℂ   — face connection-Laplacian Fiedler eigenmode (TessMat.nxr)
%
%   arg(Z_joint(x,t)) = arg(s_face(x,t)) − arg(u₁(x))   [wave phase in Fiedler gauge]
%
%   IMPORTANT — why the FULL field, not a single mode coefficient:
%     A traveling wave is a superposition of multiple spatial modes with
%     phase-shifted temporal coefficients.  Using only the dominant mode
%     coefficient θ₁(t) would make Z_joint rank-1 (= −φ_F(x) trivially) and
%     reveal no wave.  The full s_face(x,t) carries the multi-mode spatial
%     structure; the unit Fiedler phasor only fixes the spatial gauge/sign.
%
%   Sign ambiguity is removed: conj(u₁/|u₁|) absorbs the sulcal-wall π-jumps
%   carried by u₁'s gauge, not by the physical wave.
%
% INPUTS:
%   s_face      [nLHF × nTime] complex — full face source field at one scale
%   FaceIndices [nLHF × 1]    global face indices into TessMat.Faces
%   SurfaceFile Brainstorm cortex surface file path
%
% OPTIONS (name-value):
%   'ModeIndex'  which connection eigenmode is the Fiedler coordinate (default: 1)
%   'Verbose'    logical (default: true)
%
% OUTPUTS:
%   Z_joint    [nLHF × nTime] complex — gauge-demodulated face field
%   phi_F      [nLHF × 1]  real — Fiedler longitude = arg(u₁(x)) ∈ (−π,π]
%   phi_T      [1 × nTime] real — mean phase of Z_joint across active faces
%   confidence [nLHF × nTime] real — |u₁(x)| · |s_face(x,t)| product
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

nLHF  = size(s_face, 1);
nTime = size(s_face, 2);
FaceIndices = FaceIndices(:);

%% ── Load face Fiedler eigenvector from TessMat.nxr ───────────────────────
TessMat = in_tess_bst(SurfaceFile, 0);
if ~isfield(TessMat,'nxr') || isempty(TessMat.nxr)
    error('bst_fiedler_joint_field:NoNxr', ...
        'TessMat.nxr missing. Run tess_nxr_populate first.');
end
NxrData = TessMat.nxr;

% Map FaceIndices to LH-local index in ConnEigLH
lH_f_nxr = NxrData.lH_f;
[~, lhf_loc] = ismember(FaceIndices, lH_f_nxr);
valid_f = lhf_loc > 0;

u1 = zeros(nLHF, 1, 'like', 1i);
u1(valid_f) = NxrData.ConnEigLH.Vectors(lhf_loc(valid_f), ModeIndex);

phi_F = angle(u1);          % [nLHF × 1] Fiedler longitude
amp_F = abs(u1);            % confidence per face

% Unit Fiedler phasor (gauge only — amplitude removed)
u1_unit = u1 ./ max(abs(u1), eps);

%% ── Demodulate full field into Fiedler gauge ─────────────────────────────
% Z_joint(x,t) = conj(û₁(x)) · s_face(x,t)
% This rotates each face's complex value by −arg(u₁(x)), removing the spatial
% gauge (and its sulcal sign flips) while keeping the full temporal+spatial
% structure of the wave intact.
Z_joint    = conj(u1_unit) .* s_face;          % [nLHF × nTime] complex
confidence = amp_F .* abs(s_face);             % [nLHF × nTime] joint confidence

% Temporal phase: amplitude-weighted mean phase across faces at each time
phi_T = zeros(1, nTime);
w = amp_F .* (amp_F > 0);
for t = 1:nTime
    z_t = sum(w .* Z_joint(:,t)) / max(sum(w), eps);
    phi_T(t) = angle(z_t);
end

if Verbose
    fprintf('bst_fiedler_joint_field: %d LH faces  Fiedler mode=%d\n', nLHF, ModeIndex);
    fprintf('  Fiedler |u₁| range: [%.3f %.3f]  phi_F range: [%.2f %.2f] rad\n', ...
        min(amp_F), max(amp_F), min(phi_F), max(phi_F));
    fprintf('  |s_face| mean=%.3e  Z_joint demodulated (full multi-mode field)\n', ...
        mean(abs(s_face(:))));
end
end
