function WF = bst_wavefront_track(s_corr, SurfaceFile, varargin)
% BST_WAVEFRONT_TRACK  Traveling wavefront extraction from a complex source field.
%
% USAGE:
%   WF = bst_wavefront_track(s_corr, SurfaceFile)
%   WF = bst_wavefront_track(s_corr, SurfaceFile, 'Time', t, 'CenterFreq', 10, ...)
%
% DESCRIPTION:
%   Given a sign-corrected complex source field s_corr [nV x nTime] (analytic
%   signal; abs = envelope, angle = instantaneous phase), computes per-vertex
%   wave speed and direction using the surface phase gradient, plus summary
%   statistics (dominant direction, PLV, isoline segments) at every time step.
%
%   Phase gradient is estimated at each left-hemisphere vertex by projecting
%   edge-wise phase differences onto the tangent frame {e1,e2} (trivial-connection
%   frame from tess_tangents + bst_tangent_face2vertex) and solving the 2×2
%   normal equations once per vertex.  All geometry is precomputed; the time
%   loop is fully vectorised.
%
% INPUTS:
%   s_corr       [nV x nTime] complex source field (any complex source estimate)
%   SurfaceFile  Brainstorm cortex surface file path
%
% OPTIONS (name-value):
%   'Time'         [1 x nTime] time axis in seconds   (default: 1:nTime)
%   'LHOnly'       logical, restrict to left hemisphere (default: true)
%   'AmpThreshold' fraction of max amplitude for active mask (default: 0.08)
%   'CenterFreq'   Hz, used for speed = 2*pi*f0 / |grad_phi| (default: 10)
%   'IsoPhase'     rad, isoline extraction at this phase value (default: 0)
%   'Verbose'      logical, print progress (default: true)
%
% OUTPUT  WF struct:
%   .GradMag    [nV x nTime]  rad/m, NaN outside active mask
%   .Speed      [nV x nTime]  m/s,   NaN outside active mask
%   .Direction  [nV x nTime]  rad,   wave propagation direction in tangent frame
%   .DomDir     [1 x nTime]   rad,   amplitude-weighted dominant direction
%   .MeanSpeed  [1 x nTime]   m/s,   amplitude-weighted median speed
%   .PLV        [1 x nTime]   [0,1], directional phase-locking value
%   .Isolines   {nTime x 1}   cell, each [nSeg x 2 x 3] segment endpoints [m]
%   .ActiveMask [nV x 1]      logical active vertex mask
%   .e1         [nV x 3]      tangent frame e1
%   .e2         [nV x 3]      tangent frame e2
%   .Time       [1 x nTime]
%   .CenterFreq scalar Hz
%   .SurfaceFile string
%   .ComputeTime scalar s
%
% SEE ALSO: tess_tangents, bst_tangent_face2vertex, tess_hemisplit
%
% Authors: Diellor Basha, 2026

t0 = tic;

%% ── Parse options ─────────────────────────────────────────────────────────
[nV, nTime] = size(s_corr);
opts.Time         = 1:nTime;
opts.LHOnly       = true;
opts.AmpThreshold     = 0.08;   % spatial: fraction of max vertex-mean amplitude
opts.AmpThresholdTime = 0.20;   % temporal: fraction of max mean-LH amplitude
                                 % time points below this are excluded from PLV/speed/dir
opts.CenterFreq   = 10;
opts.IsoPhase     = 0;
opts.SkipIsolines = true;
opts.Verbose      = true;
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end
f0   = opts.CenterFreq;
tVec = opts.Time(:)';
if numel(tVec) ~= nTime
    error('Length of ''Time'' (%d) must match nTime (%d).', numel(tVec), nTime);
end

%% ── Load surface ──────────────────────────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile);
Vtx     = TessMat.Vertices;           % [nV x 3] metres
Faces   = double(TessMat.Faces);      % [nF x 3]

%% ── Tangent frame (trivial-connection) ────────────────────────────────────
[Uf, ~] = tess_tangents(SurfaceFile, 'NoSave', 1);
[e1, e2] = bst_tangent_face2vertex(Faces, Uf, TessMat.VertNormals);
% e1,e2: [nV x 3] per-vertex smooth tangent frame

%% ── Hemisphere mask ───────────────────────────────────────────────────────
if opts.LHOnly
    [~, lhIdx] = tess_hemisplit(TessMat);
    lhIdx = lhIdx(:);
else
    lhIdx = (1:nV)';
end
nLH = numel(lhIdx);

%% ── Amplitude mask (time-averaged envelope) ───────────────────────────────
ampEnv    = abs(s_corr);                       % [nV x nTime]
ampMean   = mean(ampEnv, 2);                   % [nV x 1]
thresh    = opts.AmpThreshold * max(ampMean(lhIdx));
activeFull = false(nV, 1);
activeFull(lhIdx) = ampMean(lhIdx) >= thresh;

%% ── Vertex connectivity (from faces) ──────────────────────────────────────
% Build sparse adjacency over lhIdx vertices only.
% Re-index to 1..nLH for geometry precomputation.
lhSet  = false(nV, 1);  lhSet(lhIdx) = true;
lhMap  = zeros(nV, 1);  lhMap(lhIdx) = 1:nLH;

% Faces fully within LH
faceInLH = lhSet(Faces(:,1)) & lhSet(Faces(:,2)) & lhSet(Faces(:,3));
F_lh     = lhMap(Faces(faceInLH, :));   % [nF_lh x 3], indices in 1..nLH
Vtx_lh   = Vtx(lhIdx, :);               % [nLH x 3]

% Neighbour lists via accumarray-style construction
adj = sparse([F_lh(:,1); F_lh(:,2); F_lh(:,3); ...
              F_lh(:,1); F_lh(:,2); F_lh(:,3)], ...
             [F_lh(:,2); F_lh(:,3); F_lh(:,1); ...
              F_lh(:,3); F_lh(:,1); F_lh(:,2)], ...
             ones(6*size(F_lh,1), 1), nLH, nLH);
adj = adj > 0;

%% ── Precompute gradient operator geometry ─────────────────────────────────
% For each LH vertex v: gather neighbours w; build A_v [k x 2] of projected
% edge vectors; form AtA_inv_v [2 x 2] and At_v [2 x k].
% Store in padded arrays — k_max neighbours.
if opts.Verbose
    fprintf('bst_wavefront_track: precomputing gradient geometry for %d LH vertices...\n', nLH);
end

e1_lh = e1(lhIdx, :);   % [nLH x 3]
e2_lh = e2(lhIdx, :);

% Max neighbour count
nbrs_all = cell(nLH, 1);
for v = 1:nLH
    nbrs_all{v} = find(adj(v,:));
end
kMax = max(cellfun(@numel, nbrs_all));

% Padded storage: At_all [nLH x 2 x kMax], AtAinv [nLH x 2 x 2]
% NaN-padded so unused entries don't contribute when summed.
At_all    = zeros(nLH, 2, kMax);   % precomputed A'  (2 x k per vertex)
AtAinv    = zeros(nLH, 2, 2);
nbrIdx    = zeros(nLH, kMax, 'uint32');   % neighbour LH-index per slot
kCount    = zeros(nLH, 1, 'uint8');

for v = 1:nLH
    nbrs = nbrs_all{v};
    k    = numel(nbrs);
    if k < 2, continue; end   % degenerate vertex: skip gradient
    kCount(v) = k;
    nbrIdx(v, 1:k) = uint32(nbrs);

    % Edge vectors projected onto tangent frame
    dXYZ = Vtx_lh(nbrs,:) - Vtx_lh(v,:);   % [k x 3]
    du   = dXYZ * e1_lh(v,:)';              % [k x 1]
    dv   = dXYZ * e2_lh(v,:)';
    A_v  = [du, dv];                         % [k x 2]

    AtA_v          = A_v' * A_v;            % [2 x 2]
    if rcond(AtA_v) < 1e-12, continue; end  % numerically degenerate
    AtAinv_v       = inv(AtA_v);            % [2 x 2]
    At_v           = A_v';                  % [2 x k]

    At_all(v, :, 1:k)  = At_v;
    AtAinv(v, :, :)    = AtAinv_v;
end

%% ── Phase-field and amplitude on LH ──────────────────────────────────────
phi_lh = angle(s_corr(lhIdx, :));   % [nLH x nTime]
amp_lh = ampEnv(lhIdx, :);          % [nLH x nTime]
activeL = activeFull(lhIdx);        % [nLH x 1]

%% ── Gradient computation (vectorised over time) ───────────────────────────
% For each vertex v and time t:
%   dphi_v [k x 1] = wrap(phi_lh(nbrs,t) - phi_lh(v,t))
%   grad_uv [2 x 1] = AtAinv_v * At_v * dphi_v
%
% We loop over time (nTime dimension) but keep nLH vectorised.

if opts.Verbose
    fprintf('bst_wavefront_track: computing phase gradient over %d time steps...\n', nTime);
end

gradU  = zeros(nLH, nTime);
gradV  = zeros(nLH, nTime);

% Precompute B_all = AtAinv * At: [nLH x 2 x kMax] — full gradient operator.
B_all = zeros(nLH, 2, kMax);
for v = 1:nLH
    k = kCount(v);
    if k < 2, continue; end
    Ai = squeeze(AtAinv(v,:,:));          % [2 x 2]
    At = squeeze(At_all(v,:,1:k));        % [2 x k]
    B_all(v,:,1:k) = Ai * At;            % [2 x k]
end

% Fix nbrIdx padding: replace 0 (unused slots) with 1 so gather never
% indexes out of bounds — B_all is 0 at those slots so contribution = 0.
nbrIdxFix = nbrIdx;
nbrIdxFix(nbrIdxFix == 0) = 1;

% Time loop — vectorised over all LH vertices via 3-D broadcast.
% phi_nbrs [nLH x kMax]: gathered neighbour phases
% dphi     [nLH x kMax]: wrapped phase differences
% prod     [nLH x 2 x kMax]: B_all .* dphi (broadcast)
% g        [nLH x 2]: summed gradient components
for t = 1:nTime
    phi_t    = phi_lh(:, t);                          % [nLH x 1]
    phi_nbrs = reshape(phi_t(nbrIdxFix(:)), nLH, kMax); % [nLH x kMax]
    dphi     = phi_nbrs - phi_t;                      % [nLH x kMax]
    dphi     = dphi - 2*pi * round(dphi / (2*pi));   % wrap to [-pi, pi]
    % B_all [nLH x 2 x kMax] times dphi [nLH x 1 x kMax] → sum over k
    g        = sum(B_all .* reshape(dphi, nLH, 1, kMax), 3); % [nLH x 2]
    gradU(:, t) = g(:, 1);
    gradV(:, t) = g(:, 2);
end

%% ── Speed, direction, masking ─────────────────────────────────────────────
gradMag = sqrt(gradU.^2 + gradV.^2);          % [nLH x nTime] rad/m
speed   = 2*pi*f0 ./ max(gradMag, eps);       % [nLH x nTime] m/s
speed   = min(speed, 50);                      % clamp to physiological range
theta   = atan2(gradV, gradU);                 % [nLH x nTime] rad

% NaN outside active mask
nanMask = ~activeL;
gradMag(nanMask, :) = NaN;
speed(nanMask, :)   = NaN;
theta(nanMask, :)   = NaN;

%% ── Temporal amplitude gate ───────────────────────────────────────────────
% Exclude time steps where mean LH amplitude is below AmpThresholdTime×peak.
% Low-amplitude periods have unreliable instantaneous phase (SNR << 1),
% corrupting PLV and wave-direction estimates.
amp_mean_t  = mean(amp_lh, 1);                         % [1 x nTime]
amp_time_thr = opts.AmpThresholdTime * max(amp_mean_t);
validTime    = amp_mean_t >= amp_time_thr;             % [1 x nTime] logical

if opts.Verbose
    fprintf('  Temporal gate (%.0f%% amp thr): %d / %d time steps active (%.1f%%)\n', ...
        100*opts.AmpThresholdTime, sum(validTime), nTime, 100*mean(validTime));
end

%% ── Dominant direction + PLV (per time step) ──────────────────────────────
domDir   = nan(1, nTime);
plv      = nan(1, nTime);
meanSpd  = nan(1, nTime);

for t = 1:nTime
    if ~validTime(t)
        continue;   % skip low-amplitude time steps
    end
    act  = activeL & kCount >= 2;
    if ~any(act)
        continue;
    end
    A_t  = amp_lh(act, t);
    th_t = theta(act, t);
    spd_t = speed(act, t);
    wSum = sum(A_t);
    if wSum < eps
        domDir(t) = NaN;  plv(t) = NaN;  meanSpd(t) = NaN;
        continue;
    end
    r          = sum(A_t .* exp(1i * th_t)) / wSum;
    domDir(t)  = angle(r);
    plv(t)     = abs(r);
    % amplitude-weighted median: sort and find weighted median
    [spd_s, si] = sort(spd_t);
    Aw   = A_t(si) / wSum;
    cw   = cumsum(Aw);
    iMed = find(cw >= 0.5, 1);
    meanSpd(t) = spd_s(max(iMed, 1));
end

%% ── Isoline extraction (per time step) ────────────────────────────────────
% Isolines are extracted in original vertex coordinates (metres).
Vtx_lh_m  = Vtx_lh;   % already in metres
activeFaceL = activeL(F_lh(:,1)) & activeL(F_lh(:,2)) & activeL(F_lh(:,3));
F_act      = F_lh(activeFaceL, :);   % [nFact x 3]

% Isoline extraction: O(nFaces×nTime) — skip by default for large recordings.
% Pass 'SkipIsolines', false to compute (useful for visualization).
Isolines = cell(nTime, 1);
if ~opts.SkipIsolines
    if opts.Verbose
        fprintf('bst_wavefront_track: extracting isolines...\n');
    end
    isoPhase = opts.IsoPhase;
    for t = 1:nTime
        Isolines{t} = local_phase_isoline(Vtx_lh_m, F_act, phi_lh(:,t), isoPhase);
    end
end

%% ── Pack into output struct (full nV space, NaN outside LH) ───────────────
GradMagOut = nan(nV, nTime);  GradMagOut(lhIdx, :) = gradMag;
SpeedOut   = nan(nV, nTime);  SpeedOut(lhIdx, :)   = speed;
DirOut     = nan(nV, nTime);  DirOut(lhIdx, :)     = theta;
e1Out      = zeros(nV, 3);    e1Out(lhIdx, :)      = e1_lh;
e2Out      = zeros(nV, 3);    e2Out(lhIdx, :)      = e2_lh;
ActiveOut  = false(nV, 1);    ActiveOut(lhIdx)     = activeL;

WF.GradMag     = GradMagOut;
WF.Speed       = SpeedOut;
WF.Direction   = DirOut;
WF.DomDir      = domDir;
WF.MeanSpeed   = meanSpd;
WF.PLV         = plv;
WF.Isolines    = Isolines;
WF.ActiveMask   = ActiveOut;
WF.ValidTime    = validTime;         % [1 x nTime] logical — high-amplitude time steps
WF.e1          = e1Out;
WF.e2          = e2Out;
WF.Time        = tVec;
WF.CenterFreq  = f0;
WF.SurfaceFile = SurfaceFile;
WF.ComputeTime = toc(t0);

if opts.Verbose
    fprintf('bst_wavefront_track: done in %.1f s.\n', WF.ComputeTime);
end

end   % bst_wavefront_track


%% ── Local helper: phase isoline extraction ────────────────────────────────
function segs = local_phase_isoline(Vtx, Faces, ph, isoPhase)
% Extract mesh segments where the phase field equals isoPhase.
% Uses zero-crossings of sin(ph - isoPhase) with cos(ph - isoPhase) > 0
% to avoid the +/-pi wrap artefact.
%
% Vtx   [nLH x 3] vertex positions
% Faces [nF x 3]  face connectivity (LH-local indices)
% ph    [nLH x 1] phase field (rad)
% Returns segs [nSeg x 2 x 3] segment endpoint coordinates.

segs  = zeros(0, 2, 3);
s     = sin(ph - isoPhase);
c     = cos(ph - isoPhase);
edges = [1 2; 2 3; 3 1];

for fi = 1:size(Faces,1)
    vv  = Faces(fi,:);
    sv  = s(vv);
    cv  = c(vv);
    pts = zeros(0, 3);
    for e = 1:3
        a = edges(e,1);  b = edges(e,2);
        if sv(a) * sv(b) < 0
            % Linear interpolation of the zero crossing
            t_cross = sv(a) / (sv(a) - sv(b));
            cmid    = cv(a) + t_cross * (cv(b) - cv(a));
            if cmid > 0
                pts(end+1,:) = Vtx(vv(a),:) + t_cross * (Vtx(vv(b),:) - Vtx(vv(a),:)); %#ok<AGROW>
            end
        end
    end
    if size(pts,1) == 2
        segs(end+1,:,:) = pts; %#ok<AGROW>
    end
end
end
