function [s_corr, sigma, info] = bst_source_sign_correct(s_field, SurfaceFile, varargin)
% BST_SOURCE_SIGN_CORRECT  Fix sulcal-wall sign ambiguity via facing-pair analysis.
%
% USAGE:
%   [s_corr, sigma, info] = bst_source_sign_correct(s_field, SurfaceFile)
%   [s_corr, sigma, info] = bst_source_sign_correct(s_field, SurfaceFile, 'MaxDist', 0.003)
%
% DESCRIPTION:
%   At sulcal walls, the cortical normal n̂ flips between facing vertices.
%   This creates π-jumps in the phase map of J·n̂ even though both walls fire
%   together — it is a sign convention artefact, not a physiological discontinuity.
%
%   Algorithm — Sulcal-pair-aware sign correction:
%
%   1. Detect facing sulcal-wall pairs (v1, v2) using sa_sulcal_walls: pairs
%      where both vertices are sulcal, geometrically close (≤MaxDist), and have
%      anti-parallel normals (n̂₁·n̂₂ < NormalDot), but are not mesh-adjacent.
%
%   2. For each pair where the phase difference |arg(ŝ(v1)/ŝ(v2))| > π/2 at
%      peak time, one vertex is sign-flipped. The vertex to flip is chosen by
%      checking GYRAL-NEIGHBOR CONTEXT: whichever wall vertex disagrees MORE
%      with its non-sulcal (gyral) neighbors is the inconsistent one and gets
%      σ=-1.  Fiedler amplitude breaks ties: flip the vertex with lower |u₁(v)|
%      (less well-defined direction → less confident measurement).
%
%   3. Aggregate votes across all pairs per vertex (majority vote weighted by
%      Fiedler amplitude of the partner vertex).
%
%   4. Optional post-pass: one round of spatial coherence to resolve isolated
%      conflicts created by step 3.
%
% WHY SPATIAL VOTING ALONE FAILS:
%   Both walls of a sulcus start at σ=+1 and their neighborhoods are
%   predominantly σ=+1, so iterated voting from a neutral state never detects
%   a problem — the error is symmetric.  Explicit pair detection breaks the
%   symmetry by comparing ACROSS the sulcus.
%
% INPUTS:
%   s_field      [nV x nTime]  complex source field (e.g. ImageGridAmp from
%                              bst_eigenmode_analytic_inverse)
%   SurfaceFile  Brainstorm relative path to the cortical surface
%
% OPTIONS (name-value):
%   'MaxDist'     0.003  m   max 3-D distance for sulcal-pair detection
%   'NormalDot'  -0.70       n̂·n̂ threshold (more negative = stricter anti-alignment)
%   'Nring'       3          mesh-adjacency ring to exclude same-wall pairs
%   'LHOnly'      true        restrict to left hemisphere
%   'SmoothPass'  false       one spatial coherence pass after pair correction
%                             (default false: smooth pass tends to revert good corrections)
%   'Verbose'     true
%
% OUTPUTS:
%   s_corr   [nV x nTime]  corrected source field
%   sigma    [nV x 1]      ±1 sign map (LH vertices only; RH = +1)
%   info     struct:
%              .nPairs       — detected sulcal wall pairs
%              .nFlipped     — vertices set to σ=-1
%              .pairPhaseDisc — fraction of pairs with |Δφ| > π/2 before correction
%              .lH_v         — LH vertex indices
%
% SEE ALSO: sa_sulcal_walls, bst_conn_eigenmodes_ensure, bst_wavefront_track
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
MaxDist    = 0.003;
NormalDot  = -0.70;
Nring      = 3;
LHOnly     = true;
SmoothPass = false;
Verbose    = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'maxdist',    MaxDist    = varargin{k+1};
        case 'normaldot',  NormalDot  = varargin{k+1};
        case 'nring',      Nring      = varargin{k+1};
        case 'lhonly',     LHOnly     = varargin{k+1};
        case 'smoothpass', SmoothPass = varargin{k+1};
        case 'verbose',    Verbose    = varargin{k+1};
    end
end

nV    = size(s_field, 1);
sigma = ones(nV, 1);

%% ── Load surface ─────────────────────────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile, 0);
Vtx  = TessMat.Vertices;
Nv   = TessMat.VertNormals;
if isfield(TessMat, 'VertConn') && ~isempty(TessMat.VertConn)
    VertConn = TessMat.VertConn;
else
    VertConn = tess_vertconn(Vtx, TessMat.Faces);
end
if isfield(TessMat, 'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = double(TessMat.SulciMap(:));
else
    SulciMap = zeros(nV, 1);
    if Verbose
        warning('bst_source_sign_correct:NoSulciMap', ...
            'SulciMap absent; pair detection disabled. Returning sigma=ones.');
    end
    s_corr = s_field;
    info   = struct('nPairs',0,'nFlipped',0,'pairPhaseDisc',NaN,'lH_v',[]);
    return;
end

%% ── LH indices ───────────────────────────────────────────────────────────
[~, lH_v] = tess_hemisplit(TessMat);
lH_v = lH_v(:);
if ~LHOnly
    lH_v = (1:nV)';
end
lhSet = false(nV,1);  lhSet(lH_v) = true;

%% ── Fiedler amplitude (confidence weights) ───────────────────────────────
ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile);
amp_F   = zeros(nV, 1);
if ~isempty(ConnEig) && isfield(ConnEig,'Vectors') && ~isempty(ConnEig.Vectors)
    if isfield(ConnEig,'Component') && isfield(ConnEig,'CompRank') && ...
            ~isempty(ConnEig.Component) && ~isempty(ConnEig.CompRank)
        fc = find(ConnEig.Component == 1 & ConnEig.CompRank == 1, 1);
    else
        fc = 1;
    end
    if isempty(fc), fc = 1; end
    zF    = ConnEig.Vectors(:, fc);
    amp_F = abs(zF);
    med_lh = median(amp_F(lH_v));
    if med_lh > eps, amp_F = amp_F / med_lh; end
end

%% ── Peak time snapshot ───────────────────────────────────────────────────
mean_amp = mean(abs(s_field(lH_v, :)), 1);
[~, t_peak] = max(mean_amp);
z_peak = s_field(:, t_peak);   % [nV x 1] complex

%% ── Detect sulcal-wall pairs (LH only) ───────────────────────────────────
pair_opts = struct('MaxDist', MaxDist, 'NormalDot', NormalDot, 'Nring', Nring);
lh_sulci  = SulciMap .* lhSet;   % restrict to LH sulcal vertices
pairs_all = tess_sulcal_pairs(Vtx, Nv, lh_sulci, VertConn, pair_opts);
nPairs    = size(pairs_all, 1);

if Verbose
    fprintf('bst_source_sign_correct: %d sulcal-wall pairs  (LH, dist≤%.1fmm)\n', ...
        nPairs, MaxDist*1000);
end

if nPairs == 0
    s_corr = s_field;
    info   = struct('nPairs',0,'nFlipped',0,'pairPhaseDisc',NaN,'lH_v',lH_v);
    if Verbose, fprintf('  No pairs — returning sigma=ones.\n'); end
    return;
end

%% ── Per-pair phase comparison and vote accumulation ──────────────────────
v1 = pairs_all(:,1);
v2 = pairs_all(:,2);

% Phase difference at peak time (wrapped to [0, π])
dphi = abs(angle(z_peak(v1) .* conj(z_peak(v2))));  % [nPairs x 1] in [0, π]
disc = dphi > pi/2;   % pairs with phase discontinuity

pairPhaseDisc = mean(disc);
if Verbose
    fprintf('  Pairs with |Δφ| > π/2: %d / %d  (%.1f%%)\n', ...
        sum(disc), nPairs, 100*pairPhaseDisc);
end

% For discontinuous pairs: determine which vertex to flip.
% vote_neg(v): weighted votes for σ(v)=-1
% vote_pos(v): weighted votes for σ(v)=+1
vote_neg = zeros(nV, 1);
vote_pos = zeros(nV, 1);

disc_v1 = v1(disc);
disc_v2 = v2(disc);

% Precompute gyral-neighbor mean phasor for each LH vertex.
% Gyral context: mean of non-sulcal mesh neighbors.
gyral_mask = (SulciMap == 0) & lhSet;   % [nV x 1]
VC_lh = double(VertConn > 0);
gyral_count = VC_lh * double(gyral_mask);      % [nV x 1] number of gyral neighbors
gyral_re    = VC_lh * (double(gyral_mask) .* real(z_peak));  % weighted sum of gyral
gyral_im    = VC_lh * (double(gyral_mask) .* imag(z_peak));
has_gyral   = gyral_count > 0;
z_gyral     = complex(gyral_re ./ max(gyral_count, eps), ...
                      gyral_im ./ max(gyral_count, eps));  % mean gyral phasor

for pi_i = 1:numel(disc_v1)
    va = disc_v1(pi_i);
    vb = disc_v2(pi_i);

    % Agreement with gyral context: Re(conj(z_gyral) * z_peak) > 0 = agrees
    if has_gyral(va) && has_gyral(vb)
        agree_a = real(conj(z_gyral(va)) * z_peak(va));
        agree_b = real(conj(z_gyral(vb)) * z_peak(vb));
        if agree_a < 0 && agree_b >= 0
            flip_v = va;  keep_v = vb;
        elseif agree_b < 0 && agree_a >= 0
            flip_v = vb;  keep_v = va;
        else
            % Tiebreaker: flip lower Fiedler amplitude
            if amp_F(va) <= amp_F(vb)
                flip_v = va;  keep_v = vb;
            else
                flip_v = vb;  keep_v = va;
            end
        end
    elseif has_gyral(va)
        agree_a = real(conj(z_gyral(va)) * z_peak(va));
        if agree_a < 0
            flip_v = va;  keep_v = vb;
        else
            flip_v = vb;  keep_v = va;
        end
    elseif has_gyral(vb)
        agree_b = real(conj(z_gyral(vb)) * z_peak(vb));
        if agree_b < 0
            flip_v = vb;  keep_v = va;
        else
            flip_v = va;  keep_v = vb;
        end
    else
        % No gyral context: Fiedler tiebreaker
        if amp_F(va) <= amp_F(vb)
            flip_v = va;  keep_v = vb;
        else
            flip_v = vb;  keep_v = va;
        end
    end

    w = amp_F(keep_v) + eps;   % weight by confidence of the anchor vertex
    vote_neg(flip_v) = vote_neg(flip_v) + w;
    vote_pos(keep_v) = vote_pos(keep_v) + w;
end

% Majority vote per vertex
sigma(lH_v) = 1;
flip_lh = lH_v(vote_neg(lH_v) > vote_pos(lH_v));
sigma(flip_lh) = -1;

if Verbose
    fprintf('  After pair votes: σ=-1 on %d LH vertices\n', numel(flip_lh));
end

%% ── Optional spatial coherence pass ─────────────────────────────────────
% One round of amplitude-weighted majority vote to smooth isolated conflicts
% caused by inconsistent pair votes on the same vertex.
if SmoothPass && numel(flip_lh) > 0
    amp_lh   = amp_F(lH_v);
    C_lh     = VertConn(lH_v, lH_v);
    sigma_lh = sigma(lH_v);
    z_lh     = z_peak(lH_v);

    weighted  = amp_lh .* sigma_lh .* z_lh;
    W_ref     = C_lh * weighted;
    W_nrm     = C_lh * amp_lh;
    z_ref     = W_ref ./ max(W_nrm, eps);
    agree_pos = real(conj(z_ref) .* z_lh);

    sigma_smooth = ones(numel(lH_v), 1);
    sigma_smooth(agree_pos < 0) = -1;
    n_change = sum(sigma_lh ~= sigma_smooth);
    sigma(lH_v) = sigma_smooth;

    if Verbose
        fprintf('  Smooth pass: changed %d vertices.  Final σ=-1: %d\n', ...
            n_change, sum(sigma(lH_v) == -1));
    end
end

%% ── Apply and pack ───────────────────────────────────────────────────────
s_corr = sigma .* s_field;

info.nPairs        = nPairs;
info.nFlipped      = sum(sigma(lH_v) == -1);
info.pairPhaseDisc = pairPhaseDisc;
info.lH_v          = lH_v;

if Verbose
    fprintf('  Done: %d LH vertices sign-corrected (%.1f%%)\n', ...
        info.nFlipped, 100*info.nFlipped/numel(lH_v));
end
end
