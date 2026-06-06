function WF = bst_face_wavefront_track(s_face, FaceIndices, SurfaceFile, varargin)
% BST_FACE_WAVEFRONT_TRACK  Traveling wavefront extraction from a complex face source field.
%
% USAGE:
%   WF = bst_face_wavefront_track(s_face, FaceIndices, SurfaceFile)
%   WF = bst_face_wavefront_track(s_face, FaceIndices, SurfaceFile, 'Time', t, ...)
%
% DESCRIPTION:
%   Face-space counterpart of bst_wavefront_track.  Operates on the DUAL MESH —
%   face centroids as nodes, primal edges connecting adjacent face centroids —
%   to maintain physical consistency with the face-based 2-form source model.
%
%   Given s_face [nLHF x nTime] (output of bst_eigenmode_analytic_inverse with
%   FaceSpace=true), computes per-face wave speed and direction, dominant direction,
%   PLV, and a Poisson-reconstructed phase potential to locate the wave source.
%
% INPUTS:
%   s_face       [nLHF x nTime] complex face source field (A·m, primal 2-form)
%   FaceIndices  [nLHF x 1]    global face indices into TessMat.Faces
%   SurfaceFile  Brainstorm cortex surface file path
%
% OPTIONS (name-value):
%   'Time'             [1 x nTime] time axis in seconds      (default: 1:nTime)
%   'AmpThreshold'     fraction of peak face-mean amplitude  (default: 0.08)
%   'AmpThresholdTime' temporal amplitude gate               (default: 0.20)
%   'CenterFreq'       Hz for speed = 2*pi*f0 / |grad_phi|  (default: 10)
%   'IsoPhase'         rad for isoline extraction            (default: 0)
%   'SkipIsolines'     logical                               (default: true)
%   'PoissonAlpha'     Tikhonov regularization for Poisson   (default: 0)
%   'Verbose'          logical, print progress               (default: true)
%
% OUTPUT  WF struct:
%   .GradMag    [nLHF x nTime]  rad/m, NaN outside active mask
%   .Speed      [nLHF x nTime]  m/s,   NaN outside active mask
%   .Direction  [nLHF x nTime]  rad,   propagation direction in face tangent frame
%   .DomDir     [1 x nTime]     rad,   amplitude-weighted dominant direction
%   .MeanSpeed  [1 x nTime]     m/s,   amplitude-weighted median speed
%   .PLV        [1 x nTime]     [0,1], directional coherence
%   .ValidTime  [1 x nTime]     logical temporal gate
%   .PoissonU   [nLHF x 1]     phase potential (at peak-PLV valid time)
%   .SourceFace scalar          global face index of wave source (argmax PoissonU)
%   .SourcePos  [1 x 3]         source face centroid [m]
%   .FaceIndices [nLHF x 1]    input FaceIndices
%   .Centroids  [nLHF x 3]     LH face centroids [m]
%   .e1         [nLHF x 3]     face tangent frame e1
%   .e2         [nLHF x 3]     face tangent frame e2
%   .Isolines   {nTime x 1}    isoline segments (if SkipIsolines=false)
%   .ActiveMask [nLHF x 1]     logical
%   .Time       [1 x nTime]
%   .CenterFreq scalar Hz
%   .SurfaceFile string
%   .ComputeTime scalar s
%
% SEE ALSO: bst_wavefront_track, tess_tangents, bst_face_eigenmode_leadfield
%
% Authors: Diellor Basha, 2026

t0 = tic;

%% ── Parse options ─────────────────────────────────────────────────────────
[nLHF, nTime] = size(s_face);
FaceIndices = FaceIndices(:);
opts.Time             = 1:nTime;
opts.AmpThreshold     = 0.08;
opts.AmpThresholdTime = 0.20;
opts.CenterFreq       = 10;
opts.IsoPhase         = 0;
opts.SkipIsolines     = true;
opts.PoissonAlpha     = 0;
opts.Verbose          = true;
for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
end
f0   = opts.CenterFreq;
tVec = opts.Time(:)';
if numel(tVec) ~= nTime
    error('Length of ''Time'' (%d) must match nTime (%d).', numel(tVec), nTime);
end

%% ── Load surface and face centroids ──────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile);
Vtx     = TessMat.Vertices;           % [nV x 3] metres
Faces   = double(TessMat.Faces);      % [nF x 3]
nF      = size(Faces, 1);

ctr_f   = (Vtx(Faces(:,1),:) + Vtx(Faces(:,2),:) + Vtx(Faces(:,3),:)) / 3; % [nF x 3]
ctr_lh  = ctr_f(FaceIndices, :);     % [nLHF x 3]

%% ── Load all operators from TessMat.nxr ─────────────────────────────────
NxrData = local_load_nxr(SurfaceFile);

e1_lh = NxrData.FaceFrames.U(FaceIndices, :);   % [nLHF × 3]
e2_lh = NxrData.FaceFrames.V(FaceIndices, :);

% Face dual Laplacian for Poisson solve (clamped positive ★₁ weights)
h1d    = abs(full(diag(NxrData.hodge1)));
nE_dec = size(NxrData.d1, 2);
h1inv  = spdiags(1./max(h1d, 1e-10*max(h1d)), 0, nE_dec, nE_dec);
L_lh   = (NxrData.d1 * h1inv * NxrData.d1');
L_lh   = L_lh(FaceIndices, FaceIndices);

FA_lh  = NxrData.FaceAdj(FaceIndices, FaceIndices);

%% ── Amplitude mask ────────────────────────────────────────────────────────
ampEnv   = abs(s_face);                             % [nLHF x nTime]
ampMean  = mean(ampEnv, 2);                         % [nLHF x 1]
thresh   = opts.AmpThreshold * max(ampMean);
activeF  = ampMean >= thresh;                       % [nLHF x 1] logical

amp_mean_t = mean(ampEnv, 1);                       % [1 x nTime]
validTime  = amp_mean_t >= opts.AmpThresholdTime * max(amp_mean_t);

if opts.Verbose
    fprintf('bst_face_wavefront_track: precomputing dual-mesh gradient geometry for %d LH faces...\n', nLHF);
end

%% ── Precompute dual-mesh gradient operator geometry ──────────────────────
nbrs_all = cell(nLHF, 1);
for f = 1:nLHF
    nbrs_all{f} = find(FA_lh(f,:));
end
kMax = max(cellfun(@numel, nbrs_all));
if kMax == 0
    error('bst_face_wavefront_track: no face adjacency found. Check FaceIndices.');
end

At_all  = zeros(nLHF, 2, kMax);
AtAinv  = zeros(nLHF, 2, 2);
nbrIdx  = zeros(nLHF, kMax, 'uint32');
kCount  = zeros(nLHF, 1, 'uint8');

for f = 1:nLHF
    nbrs = nbrs_all{f};
    k    = numel(nbrs);
    if k < 2, continue; end
    kCount(f) = k;
    nbrIdx(f, 1:k) = uint32(nbrs);

    dXYZ = ctr_lh(nbrs,:) - ctr_lh(f,:);
    A_f  = [dXYZ * e1_lh(f,:)', dXYZ * e2_lh(f,:)'];  % [k x 2]
    AtA_f = A_f' * A_f;
    if rcond(AtA_f) < 1e-12, continue; end
    AtAinv(f, :, :)   = inv(AtA_f);
    At_all(f, :, 1:k) = A_f';
end

B_all = zeros(nLHF, 2, kMax);   % B = AtAinv * At  [nLHF x 2 x kMax]
for f = 1:nLHF
    k = kCount(f);
    if k < 2, continue; end
    B_all(f,:,1:k) = squeeze(AtAinv(f,:,:)) * squeeze(At_all(f,:,1:k));
end

nbrIdxFix = nbrIdx;
nbrIdxFix(nbrIdxFix == 0) = 1;   % safe padding: B_all is 0 at unused slots

%% ── Phase gradient time loop (vectorised over faces) ─────────────────────
if opts.Verbose
    fprintf('bst_face_wavefront_track: computing phase gradient over %d time steps...\n', nTime);
end

phi_lhf = angle(s_face);              % [nLHF x nTime]
gradU   = zeros(nLHF, nTime);
gradV   = zeros(nLHF, nTime);

for t = 1:nTime
    phi_t    = phi_lhf(:, t);
    phi_nbrs = reshape(phi_t(nbrIdxFix(:)), nLHF, kMax);   % [nLHF x kMax]
    dphi     = phi_nbrs - phi_t;
    dphi     = dphi - 2*pi * round(dphi / (2*pi));         % wrap to [-pi, pi]
    g        = sum(B_all .* reshape(dphi, nLHF, 1, kMax), 3); % [nLHF x 2]
    gradU(:, t) = g(:, 1);
    gradV(:, t) = g(:, 2);
end

%% ── Speed, direction, masking ─────────────────────────────────────────────
gradMag = sqrt(gradU.^2 + gradV.^2);
speed   = 2*pi*f0 ./ max(gradMag, eps);
speed   = min(speed, 50);
theta   = atan2(gradV, gradU);

nanMask = ~activeF;
gradMag(nanMask, :) = NaN;
speed(nanMask, :)   = NaN;
theta(nanMask, :)   = NaN;

if opts.Verbose
    fprintf('  Temporal gate (%.0f%% amp thr): %d / %d time steps active (%.1f%%)\n', ...
        100*opts.AmpThresholdTime, sum(validTime), nTime, 100*mean(validTime));
end

%% ── Dominant direction + PLV (per time step) ──────────────────────────────
domDir  = nan(1, nTime);
plv     = nan(1, nTime);
meanSpd = nan(1, nTime);

for t = 1:nTime
    if ~validTime(t), continue; end
    act = activeF & kCount >= 2;
    if ~any(act), continue; end
    A_t   = ampEnv(act, t);
    th_t  = theta(act, t);
    spd_t = speed(act, t);
    wSum  = sum(A_t);
    if wSum < eps
        continue;
    end
    r         = sum(A_t .* exp(1i * th_t)) / wSum;
    domDir(t) = angle(r);
    plv(t)    = abs(r);
    [spd_s, si] = sort(spd_t);
    Aw   = A_t(si) / wSum;
    cw   = cumsum(Aw);
    iMed = find(cw >= 0.5, 1);
    meanSpd(t) = spd_s(max(iMed, 1));
end

%% ── Poisson source solve (at peak-PLV valid time) ─────────────────────────
plv_valid = plv;  plv_valid(~validTime) = -Inf;
[~, t_poisson] = max(plv_valid);
if isinf(plv_valid(t_poisson))
    t_poisson = find(validTime, 1);
end

% LH primal edges (upper-triangle of FA_lh).
[rowE, colE] = find(triu(FA_lh));
nE_lh = numel(rowE);

% Build 1-form: project mean face gradient onto each dual-edge direction.
F_1form = zeros(nE_lh, 1);
for ei = 1:nE_lh
    f1 = rowE(ei);  f2 = colE(ei);
    d12 = ctr_lh(f2,:) - ctr_lh(f1,:);
    dn  = norm(d12);
    if dn < eps, continue; end
    g1 = gradU(f1,t_poisson)*e1_lh(f1,:) + gradV(f1,t_poisson)*e2_lh(f1,:);
    g2 = gradU(f2,t_poisson)*e1_lh(f2,:) + gradV(f2,t_poisson)*e2_lh(f2,:);
    F_1form(ei) = dot((g1+g2)/2, d12/dn) * dn;
end

% Poisson solve: (L_conn_lh + reg·I) u = d1_lh · H1inv · F_1form
% Use face connection Laplacian from TessMat.nxr (more principled than dual scalar Laplacian)
L_conn_lh = NxrData.ConnLaplacian(FaceIndices, FaceIndices);

d1_lh    = sparse([rowE; colE], repmat((1:nE_lh)',2,1), ...
                  [-ones(nE_lh,1); ones(nE_lh,1)], nLHF, nE_lh);
h1_lh    = max(sqrt(sum((ctr_lh(colE,:)-ctr_lh(rowE,:)).^2,2)), 1e-10);
H1inv_lh = spdiags(1./h1_lh, 0, nE_lh, nE_lh);

rhs = d1_lh * (H1inv_lh * F_1form);
reg = max(real(opts.PoissonAlpha), 1e-10);
u   = (L_conn_lh + reg*speye(nLHF)) \ rhs;

u_active = u;  u_active(~activeF) = -Inf;
[~, src_lhf_idx] = max(u_active);
src_face_global  = FaceIndices(src_lhf_idx);
src_pos          = ctr_f(src_face_global, :);

%% ── Isoline extraction (optional) ────────────────────────────────────────
Isolines = cell(nTime, 1);
if ~opts.SkipIsolines
    if opts.Verbose, fprintf('bst_face_wavefront_track: extracting isolines...\n'); end
    % Build centroid triangulation: triples mutually adjacent in FA_lh.
    [r3, c3] = find(tril(FA_lh));
    triList  = zeros(0, 3, 'uint32');
    for ei = 1:numel(r3)
        f1 = r3(ei);  f2 = c3(ei);
        common = find(FA_lh(f1,:) & FA_lh(f2,:));
        for fk = common
            if fk > f2
                triList(end+1,:) = uint32([f1, f2, fk]); %#ok<AGROW>
            end
        end
    end
    if ~isempty(triList)
        for t = 1:nTime
            Isolines{t} = local_phase_isoline(ctr_lh, double(triList), phi_lhf(:,t), opts.IsoPhase);
        end
    end
end

%% ── Output struct ─────────────────────────────────────────────────────────
WF.GradMag    = gradMag;
WF.Speed      = speed;
WF.Direction  = theta;
WF.DomDir     = domDir;
WF.MeanSpeed  = meanSpd;
WF.PLV        = plv;
WF.ValidTime  = validTime;
WF.PoissonU   = u;
WF.SourceFace = src_face_global;
WF.SourcePos  = src_pos;
WF.FaceIndices = FaceIndices;
WF.Centroids  = ctr_lh;
WF.e1         = e1_lh;
WF.e2         = e2_lh;
WF.Isolines   = Isolines;
WF.ActiveMask = activeF;
WF.Time       = tVec;
WF.CenterFreq = f0;
WF.SurfaceFile = SurfaceFile;
WF.ComputeTime = toc(t0);

if opts.Verbose
    fprintf('bst_face_wavefront_track: done in %.1f s.  Source face %d at [%.1f %.1f %.1f] mm.\n', ...
        WF.ComputeTime, WF.SourceFace, WF.SourcePos * 1000);
end

end   % bst_face_wavefront_track


%% ── Local helper: load or populate TessMat.nxr geometry ──────────────────
function NxrData = local_load_nxr(SurfaceFile)
    TessMat = in_tess_bst(SurfaceFile, 0);
    if isfield(TessMat,'nxr') && ~isempty(TessMat.nxr)
        NxrData = TessMat.nxr;
    else
        NxrData = tess_nxr_populate(SurfaceFile);
    end
end

%% ── Local helper: phase isoline on centroid triangulation ─────────────────
function segs = local_phase_isoline(Vtx, Faces, ph, isoPhase)
% Extract mesh segments where phase field equals isoPhase.
% Uses zero-crossings of sin(ph - isoPhase) with cos > 0 to avoid wrap artefact.
% Vtx   [nF x 3] centroid positions, Faces [nTri x 3] centroid triangle indices.
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
