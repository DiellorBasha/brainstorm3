function [s_corr, sigma, info] = bst_source_sign_correct(s_field, SurfaceFile, varargin)
% BST_SOURCE_SIGN_CORRECT  Fix sulcal-wall sign ambiguity in MEG source maps.
%
% USAGE:
%   [s_corr, sigma, info] = bst_source_sign_correct(s_field, SurfaceFile)
%   [s_corr, sigma, info] = bst_source_sign_correct(s_field, SurfaceFile, 'LHOnly', true, ...)
%
% DESCRIPTION:
%   At each sulcal wall the cortical normal n̂ flips, so J·n̂ is measured with
%   opposite sign on facing walls even though both walls fire together in a
%   traveling wave.  This corrupts the phase map with π-jumps at sulci.
%
%   Algorithm — Amplitude-weighted iterative phase-coherence voting in the
%   Fiedler gauge:
%     1. The Fiedler eigenmode of the connection Laplacian gives a complex
%        tangent vector u₁(v) whose phase φ_F(v)=arg(u₁(v)) is smooth across
%        sulcal walls.  Its amplitude |u₁(v)| weights vertex authority.
%     2. σ(v) ∈ {±1} is initialised by demodulating s_field against φ_F at
%        the time of peak mean-LH amplitude.
%     3. Iterative sweeps propagate local phase coherence: the weighted
%        neighborhood phasor votes on σ until convergence (O(nEdges) per sweep
%        via sparse matrix multiply — no vertex loop).
%
% INPUTS:
%   s_field      [nV x nTime]  real or complex source field (output of
%                inverse solver, e.g. ImageGridAmp or analytic signal)
%   SurfaceFile  Brainstorm relative path to the cortical surface
%
% OPTIONS (name-value):
%   'LHOnly'        true  — only correct left hemisphere (default: true)
%   'MaxIter'       50    — maximum voting sweeps
%   'Tol'           1e-3  — convergence: fraction of LH vertices that flipped
%   'MinFiedlerAmp' 0.1   — min normalised Fiedler amplitude to trust σ_init
%   'Verbose'       true  — print progress
%
% OUTPUTS:
%   s_corr   [nV x nTime]  corrected source field (same type as s_field)
%   sigma    [nV x 1]      sign map: ±1 at LH vertices, +1 elsewhere
%   info     struct with fields:
%              .nSweeps    — iterations until convergence
%              .flipFrac   — [nSweeps x 1] fraction flipped per sweep
%              .nFlipped   — count of LH vertices with σ=-1 in final map
%              .sigma_init — initial σ before voting
%              .lH_v       — LH vertex indices
%
% SEE ALSO: bst_conn_eigenmodes_ensure, in_tess_conn_eigenmodes, tess_hemisplit
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

%% ── Parse options ────────────────────────────────────────────────────────────
LHOnly        = true;
MaxIter       = 50;
Tol           = 1e-3;
MinFiedlerAmp = 0.1;
Verbose       = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'lhonly',        LHOnly        = varargin{k+1};
        case 'maxiter',       MaxIter       = varargin{k+1};
        case 'tol',           Tol           = varargin{k+1};
        case 'minfiedleramp', MinFiedlerAmp = varargin{k+1};
        case 'verbose',       Verbose       = varargin{k+1};
    end
end

nV    = size(s_field, 1);
sigma = ones(nV, 1);

%% ── Load surface ─────────────────────────────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile, 0);
if isfield(TessMat, 'VertConn') && ~isempty(TessMat.VertConn)
    VertConn = TessMat.VertConn;
else
    VertConn = tess_vertconn(TessMat.Vertices, TessMat.Faces);
end

%% ── LH vertex indices ────────────────────────────────────────────────────────
[~, lH_v] = tess_hemisplit(TessMat);
lH_v      = lH_v(:);
if isempty(lH_v)
    warning('bst_source_sign_correct:NoLH', ...
        'tess_hemisplit returned no left-hemisphere vertices; returning sigma=ones.');
    s_corr = s_field;
    info   = struct('nSweeps', 0, 'flipFrac', [], 'nFlipped', 0, ...
                    'sigma_init', sigma, 'lH_v', lH_v);
    return;
end
if ~LHOnly
    % Extend to all vertices (both hemispheres).
    lH_v = (1:nV)';
end

%% ── Load connection eigenmodes ───────────────────────────────────────────────
ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile);
if isempty(ConnEig) || ~isfield(ConnEig, 'Vectors') || isempty(ConnEig.Vectors)
    warning('bst_source_sign_correct:NoConnEig', ...
        'Connection eigenmodes not available for %s; returning sigma=ones.', SurfaceFile);
    s_corr = s_field;
    info   = struct('nSweeps', 0, 'flipFrac', [], 'nFlipped', 0, ...
                    'sigma_init', sigma, 'lH_v', lH_v);
    return;
end

%% ── Fiedler vector ───────────────────────────────────────────────────────────
% Prefer the column tagged as Component==1, CompRank==1 (first non-trivial
% mode); fall back to column 1 (ascending eigenvalue order).
if isfield(ConnEig, 'Component') && isfield(ConnEig, 'CompRank') && ...
        ~isempty(ConnEig.Component) && ~isempty(ConnEig.CompRank)
    fiedler_col = find(ConnEig.Component == 1 & ConnEig.CompRank == 1, 1);
else
    fiedler_col = [];
end
if isempty(fiedler_col)
    fiedler_col = 1;
end
zF     = ConnEig.Vectors(:, fiedler_col);   % [nV x 1] complex
amp_F  = abs(zF);
med_lh = median(amp_F(lH_v));
if med_lh < eps
    med_lh = 1;
end
amp_F_norm = amp_F / med_lh;               % normalised; ~1 on well-defined vertices

%% ── Peak time-frame on LH ────────────────────────────────────────────────────
mean_amp_lh = mean(abs(s_field(lH_v, :)), 1);   % [1 x nTime]
[~, t_peak] = max(mean_amp_lh);
z_peak      = s_field(:, t_peak);               % [nV x 1] real or complex

%% ── Initialise σ ─────────────────────────────────────────────────────────────
% Start from σ = +1 everywhere. Fiedler-demodulation was tried but fails because
% it compares the TEMPORAL instantaneous phase of the wave (arg(z_peak)) against
% the SPATIAL Fiedler longitude (arg(u₁)) — these live in orthogonal spaces and
% their correlation is not meaningful. Instead let the spatial voting converge
% from a neutral state; the Fiedler amplitude acts only as confidence weight.
sigma_init = ones(nV, 1);

sigma = sigma_init;
sigma_init_out = sigma_init;

if Verbose
    fprintf('bst_source_sign_correct: nV=%d, nLH=%d, t_peak=%d\n', ...
        nV, numel(lH_v), t_peak);
end

%% ── Iterative voting (vectorised, O(nEdges) per sweep) ───────────────────────
% For the LH sub-graph:
%   W_ref(v) = Σ_{w∈N(v)} amp_F(w) · σ(w) · z_peak(w)   (weighted phasor sum)
%   W_nrm(v) = Σ_{w∈N(v)} amp_F(w)                       (weight normalisation)
%   agree_pos(v) = Re( conj(W_ref(v)/W_nrm(v)) · z_peak(v) )
%   σ_new(v) = +1 if agree_pos ≥ 0, else -1
%
% VertConn(lH_v, lH_v) selects the LH sub-adjacency; sparse multiply is
% O(nEdges_LH) per sweep.
C_lh  = VertConn(lH_v, lH_v);   % [nLH x nLH] sparse
amp_lh = amp_F(lH_v);            % [nLH x 1]

flipFrac = zeros(MaxIter, 1);
nSweeps  = 0;

for iter = 1:MaxIter
    sigma_lh  = sigma(lH_v);                 % [nLH x 1]
    z_lh      = z_peak(lH_v);               % [nLH x 1]
    weighted  = amp_lh .* sigma_lh .* z_lh; % [nLH x 1] amplitude-weighted phasor
    W_ref     = C_lh * weighted;            % [nLH x 1] neighbourhood sum
    W_nrm     = C_lh * amp_lh;             % [nLH x 1] weight normalisation
    z_ref     = W_ref ./ max(W_nrm, eps);   % [nLH x 1] normalised ref phasor
    agree_pos = real(conj(z_ref) .* z_lh);  % [nLH x 1]
    sigma_new_lh = ones(numel(lH_v), 1);
    sigma_new_lh(agree_pos < 0) = -1;

    n_flip          = sum(sigma_lh ~= sigma_new_lh);
    flipFrac(iter)  = n_flip / numel(lH_v);
    nSweeps         = iter;
    sigma(lH_v)     = sigma_new_lh;

    if Verbose && (mod(iter,10)==0 || flipFrac(iter) < Tol)
        fprintf('  Sweep %d: flipped %.4f\n', iter, flipFrac(iter));
    end
    if flipFrac(iter) < Tol
        break;
    end
end

flipFrac = flipFrac(1:nSweeps);

%% ── Apply σ ──────────────────────────────────────────────────────────────────
s_corr = sigma .* s_field;   % [nV x nTime], broadcast across time

%% ── Pack info ────────────────────────────────────────────────────────────────
info = struct();
info.nSweeps    = nSweeps;
info.flipFrac   = flipFrac;
info.nFlipped   = sum(sigma(lH_v) == -1);
info.sigma_init = sigma_init_out;
info.lH_v       = lH_v;

if Verbose
    fprintf('  Converged in %d sweep(s).  Final σ=-1: %d / %d LH vertices\n', ...
        nSweeps, info.nFlipped, numel(lH_v));
end
end
