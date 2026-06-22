function [s_corr, sigma_f, info] = bst_face_sign_correct(s_face, FaceIndices, SurfaceFile, varargin)
% BST_FACE_SIGN_CORRECT  Fix sulcal-wall sign ambiguity in face-based source maps.
%
% USAGE:
%   [s_corr, sigma_f, info] = bst_face_sign_correct(s_face, FaceIndices, SurfaceFile)
%
% DESCRIPTION:
%   Face-space counterpart of bst_source_sign_correct.  Operates directly on
%   the face 2-form source field s_f = Ψ·θ̂ using the face-based connection
%   Laplacian Fiedler vector from TessMat.nxr — no vertex-to-face interpolation.
%
%   The face Fiedler vector u₁(f) ∈ ℂ is the first non-trivial eigenfunction
%   of the face connection Laplacian (from tess_nxr_populate).  Its amplitude
%   |u₁(f)| is large on gyral crown faces and near-zero at Fiedler singularities
%   (sulcal floors).  The sulcal WALL FACES — those whose three vertices all
%   lie in sulcal territory — are the candidates for sign correction.
%
%   Algorithm (same structure as bst_source_sign_correct):
%   1. Detect facing sulcal-wall FACE PAIRS using face centroid proximity,
%      anti-parallel face normals, and non-adjacency in face mesh.
%   2. For each pair with |Δφ| > π/2: choose which face to flip using gyral
%      neighbour context; Fiedler amplitude |u₁(f)| breaks ties.
%   3. Majority-vote across all pairs per face.
%
% INPUTS:
%   s_face       [nLHF × nTime]  complex face source field (A·m, primal 2-form)
%   FaceIndices  [nLHF × 1]     global face indices into TessMat.Faces
%   SurfaceFile  Brainstorm cortex surface file path
%
% OPTIONS (name-value):
%   'MaxDist'     0.003   m    max 3-D centroid distance for pair detection
%   'NormalDot'  -0.70        n̂·n̂ threshold (anti-parallel normals)
%   'Nring'       3            exclude face-adjacency within N rings
%   'SmoothPass'  false        one coherence pass after pair voting
%   'Verbose'     true
%
% OUTPUTS:
%   s_corr   [nLHF × nTime]  sign-corrected face source field
%   sigma_f  [nLHF × 1]      ±1 sign map per face
%   info     struct: .nPairs, .nFlipped, .pairPhaseDisc, .FaceIndices
%
% SEE ALSO: bst_source_sign_correct, tess_nxr_populate, bst_face_wavefront_track
%
% Authors: Diellor Basha, 2026

%% ── Parse options ────────────────────────────────────────────────────────
MaxDist    = 0.003;
NormalDot  = -0.70;
Nring      = 3;
SmoothPass = false;
Verbose    = true;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'maxdist',    MaxDist    = varargin{k+1};
        case 'normaldot',  NormalDot  = varargin{k+1};
        case 'nring',      Nring      = varargin{k+1};
        case 'smoothpass', SmoothPass = logical(varargin{k+1});
        case 'verbose',    Verbose    = varargin{k+1};
    end
end

nLHF = size(s_face, 1);
FaceIndices = FaceIndices(:);
sigma_f = ones(nLHF, 1);

%% ── Load surface + nxr data ──────────────────────────────────────────────
TessMat = in_tess_bst(SurfaceFile, 0);
Vtx     = TessMat.Vertices;
Faces   = double(TessMat.Faces);
nF_all  = size(Faces, 1);

NxrData = local_load_nxr(SurfaceFile);

% Face centroids and normals for the LH face set
ctr_all = (Vtx(Faces(:,1),:) + Vtx(Faces(:,2),:) + Vtx(Faces(:,3),:)) / 3;
ctr_lh  = ctr_all(FaceIndices, :);   % [nLHF × 3]

e1_lh = NxrData.FaceFrames.U(FaceIndices, :);
e2_lh = NxrData.FaceFrames.V(FaceIndices, :);
nhat_lh = NxrData.FaceFrames.N(FaceIndices, :);  % inward normals (nxr CW convention)

%% ── Face Fiedler vector (no interpolation — computed on faces by nxr) ────
% ConnEigLH.Vectors [nLHF_nxr × K] for all LH faces.
% Map our FaceIndices to nxr's LH face set.
lH_f_nxr = NxrData.lH_f;   % global LH face indices from nxr
[~, lhf_loc] = ismember(FaceIndices, lH_f_nxr);
valid_fiedler = lhf_loc > 0;

zF_f = zeros(nLHF, 1, 'like', 1i);
if any(valid_fiedler) && ~isempty(NxrData.ConnEigLH.Vectors)
    zF_f(valid_fiedler) = NxrData.ConnEigLH.Vectors(lhf_loc(valid_fiedler), 1);
end
amp_F_f = abs(zF_f);   % Fiedler amplitude per face: 0 near singularities
med_amp = median(amp_F_f(amp_F_f > 0));
if med_amp < eps, med_amp = 1; end
amp_F_norm = amp_F_f / med_amp;

%% ── SulciMap per face ────────────────────────────────────────────────────
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    sulci_v = double(TessMat.SulciMap(:));
    % A face is "sulcal" if majority of its vertices are sulcal
    s_sum = sulci_v(Faces(FaceIndices,1)) + sulci_v(Faces(FaceIndices,2)) + ...
            sulci_v(Faces(FaceIndices,3));
    sulci_f = double(s_sum >= 2);   % [nLHF × 1]
else
    sulci_f = ones(nLHF, 1);   % treat all faces as sulcal if map missing
end

%% ── Peak time snapshot ───────────────────────────────────────────────────
mean_amp = mean(abs(s_face), 2);
[~, t_peak] = max(mean(abs(s_face), 1));
z_peak = s_face(:, t_peak);   % [nLHF × 1] complex

%% ── Detect sulcal-wall face pairs ────────────────────────────────────────
iSulci = find(sulci_f > 0);
if isempty(iSulci)
    s_corr = s_face;
    info   = struct('nPairs',0,'nFlipped',0,'pairPhaseDisc',NaN,'FaceIndices',FaceIndices);
    if Verbose, fprintf('bst_face_sign_correct: no sulcal faces — skipping.\n'); end
    return;
end

ctr_sulci = ctr_lh(iSulci, :);
nb = rangesearch(ctr_sulci, ctr_sulci, MaxDist);

% Face adjacency (N-ring) for exclusion
FA_lh = NxrData.FaceAdj(FaceIndices, FaceIndices);
FA_ring = FA_lh;  P = FA_lh;
for r = 2:Nring
    P = double((P * FA_lh) > 0);
    FA_ring = double((FA_ring + P) > 0);
end

I_pairs = []; J_pairs = [];
for a = 1:numel(iSulci)
    gi = iSulci(a);
    for b = nb{a}(:)'
        gj = iSulci(b);
        if gj <= gi, continue; end
        if dot(nhat_lh(gi,:), nhat_lh(gj,:)) >= NormalDot, continue; end
        if FA_ring(gi, gj) ~= 0, continue; end
        I_pairs(end+1,1) = gi; J_pairs(end+1,1) = gj; %#ok<AGROW>
    end
end
nPairs = numel(I_pairs);

if Verbose
    fprintf('bst_face_sign_correct: %d sulcal face pairs  (dist≤%.1fmm)\n',...
        nPairs, MaxDist*1000);
end
if nPairs == 0
    s_corr = s_face;
    info   = struct('nPairs',0,'nFlipped',0,'pairPhaseDisc',NaN,'FaceIndices',FaceIndices);
    return;
end

%% ── Phase comparison and pair voting ─────────────────────────────────────
dphi      = abs(angle(z_peak(I_pairs) .* conj(z_peak(J_pairs))));
disc      = dphi > pi/2;
pairDisc  = mean(disc);

if Verbose
    fprintf('  Pairs |Δφ|>π/2: %d/%d (%.1f%%)\n', sum(disc), nPairs, 100*pairDisc);
end

% Gyral-neighbour context using face adjacency + Fiedler weighting
% Gyral faces = low SulciMap score (invert: sulci_f==0)
gyral_mask_f = (sulci_f == 0);
VC_f = double(FA_lh > 0);
gyral_count = VC_f * double(gyral_mask_f);
gyral_re    = VC_f * (double(gyral_mask_f) .* real(z_peak));
gyral_im    = VC_f * (double(gyral_mask_f) .* imag(z_peak));
has_gyral   = gyral_count > 0;
z_gyral     = complex(gyral_re ./ max(gyral_count, eps), ...
                      gyral_im ./ max(gyral_count, eps));

vote_neg = zeros(nLHF, 1);
vote_pos = zeros(nLHF, 1);

disc_i = I_pairs(disc);
disc_j = J_pairs(disc);

for pi_k = 1:numel(disc_i)
    fa = disc_i(pi_k); fb = disc_j(pi_k);
    if has_gyral(fa) && has_gyral(fb)
        ag_a = real(conj(z_gyral(fa)) * z_peak(fa));
        ag_b = real(conj(z_gyral(fb)) * z_peak(fb));
        if ag_a < 0 && ag_b >= 0
            flip_v = fa; keep_v = fb;
        elseif ag_b < 0 && ag_a >= 0
            flip_v = fb; keep_v = fa;
        else
            if amp_F_norm(fa) <= amp_F_norm(fb)
                flip_v = fa; keep_v = fb;
            else
                flip_v = fb; keep_v = fa;
            end
        end
    elseif has_gyral(fa)
        flip_v = fb; keep_v = fa;
    elseif has_gyral(fb)
        flip_v = fa; keep_v = fb;
    else
        if amp_F_norm(fa) <= amp_F_norm(fb)
            flip_v = fa; keep_v = fb;
        else
            flip_v = fb; keep_v = fa;
        end
    end
    w = amp_F_norm(keep_v) + eps;
    vote_neg(flip_v) = vote_neg(flip_v) + w;
    vote_pos(keep_v) = vote_pos(keep_v) + w;
end

sigma_f(vote_neg > vote_pos) = -1;

if Verbose
    fprintf('  After pair votes: σ=-1 on %d faces\n', sum(sigma_f == -1));
end

%% ── Optional spatial coherence pass ─────────────────────────────────────
if SmoothPass && any(sigma_f == -1)
    u_lh  = z_peak;
    C_lh  = FA_lh;
    w_lh  = amp_F_norm;
    u_unit = u_lh ./ max(abs(u_lh), eps);
    weighted = w_lh .* sigma_f .* u_unit;
    W_ref    = C_lh * weighted;
    W_nrm    = C_lh * w_lh;
    u_ref    = W_ref ./ max(abs(W_ref), eps);
    agree    = real(conj(u_ref) .* u_unit);
    sigma_f  = sign(agree);
    sigma_f(sigma_f == 0) = 1;
    if Verbose
        fprintf('  Smooth pass: final σ=-1 on %d faces\n', sum(sigma_f == -1));
    end
end

%% ── Apply and pack ───────────────────────────────────────────────────────
s_corr = sigma_f .* s_face;

info.nPairs        = nPairs;
info.nFlipped      = sum(sigma_f == -1);
info.pairPhaseDisc = pairDisc;
info.FaceIndices   = FaceIndices;

if Verbose
    fprintf('  Done: %d faces corrected (%.1f%%)\n',...
        info.nFlipped, 100*info.nFlipped/nLHF);
end
end

%% ── Local helper ─────────────────────────────────────────────────────────
function NxrData = local_load_nxr(SurfaceFile)
    TessMat = in_tess_bst(SurfaceFile, 0);
    if isfield(TessMat,'nxr') && ~isempty(TessMat.nxr)
        NxrData = TessMat.nxr;
    else
        NxrData = tess_nxr_populate(SurfaceFile);
    end
end
