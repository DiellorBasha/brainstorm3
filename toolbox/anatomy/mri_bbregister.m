function [TransfWorld, cost, info] = mri_bbregister(sMriPet, sMriRef, WhiteSurfFile, TinitWorld, Opts)
% MRI_BBREGISTER: Boundary-Based Registration of a (PET) volume to anatomy (Greve 2009).
%
% Native replication of FreeSurfer bbregister: finds the 6-DOF RIGID transform that aligns
% the moving volume (PET) to the subject's anatomy by MAXIMIZING the GM/WM intensity contrast
% sampled across the white surface - instead of a global intensity similarity. This is robust
% for low-contrast PET where intensity-based coregistration drifts. It is a LOCAL refiner:
% pass a coarse initial alignment and it polishes it.
%
% SURFACE-CONSTRAINED (tangential) SMOOTHING (Opts.HeatT > 0): the GM-side and WM-side
% intensity samples are fields on the white-surface vertices; they are denoised by a
% Laplace-Beltrami HEAT diffusion ALONG the surface (implicit Euler g=(M+t.K)\(M.f), cached
% Cholesky via tess_operators) - i.e. smoothing tangentially within each depth shell WITHOUT
% blurring radially across the GM/WM boundary. This improves the SNR/robustness of the cost
% landscape while preserving the radial contrast BBR relies on (unlike isotropic 3-D smoothing).
%
% PER-FRAME (4-D PET): if sMriPet.Cube is 4-D, all frames are sampled at the same pose.
% FrameCombine='mean' averages them (equivalent to smoothing the static, for linear smoothing);
% 'median' smooths each frame and takes the per-vertex median contrast - robust to motion-/
% noise-corrupted frames (a genuine per-frame gain).
%
% USAGE:  [TransfWorld, cost, info] = mri_bbregister(sMriPet, sMriRef, WhiteSurfFile, TinitWorld, Opts)
%
% INPUTS:
%   sMriPet      : moving MRI/PET struct (Cube 3-D or 4-D + SCS/world geometry).
%   sMriRef      : reference anatomical MRI struct (defines SCS<->world for the surface).
%   WhiteSurfFile: white (GM/WM) surface file, registered to sMriRef.
%   TinitWorld   : 4x4 initial PET-world -> anat-world transform ([] = identity).
%   Opts         : .SampleDist (mm, 1.5), .SubSample (max verts when no heat, 5000),
%                  .Mslope (0.5), .HeatT (heat time, 0 = off), .FrameCombine ('mean'|'median').
%
% OUTPUTS:
%   TransfWorld  : 4x4 refined PET-world -> anat-world rigid transform.
%   cost         : final BBR cost (0 = sharp boundary, 1 = no contrast).
%   info         : struct with .theta .costInit .nVert .sign .cost .heatT .nFrames.
%
% Reference: Greve DN, Fischl B. Accurate and robust brain image alignment using boundary-based
%            registration. NeuroImage 2009;48(1):63-72.
%
% SEE ALSO: mri_coregister, cs_convert, tess_operators
%
% Author: Diellor Basha, 2026

    if (nargin < 4), TinitWorld = []; end
    if isempty(TinitWorld), TinitWorld = eye(4); end
    if (nargin < 5) || isempty(Opts), Opts = struct(); end
    Def = struct('SampleDist',1.5, 'SubSample',5000, 'Mslope',0.5, 'HeatT',0, 'FrameCombine','mean');
    fn = fieldnames(Def); for i=1:numel(fn), if ~isfield(Opts,fn{i})||isempty(Opts.(fn{i})), Opts.(fn{i})=Def.(fn{i}); end; end

    % World units per millimetre (Brainstorm 'world' may be metres).
    p0 = cs_convert(sMriRef,'voxel','world',[1 1 1]); p1 = cs_convert(sMriRef,'voxel','world',[2 1 1]);
    wpm = norm(p1 - p0) / sMriRef.Voxsize(1);

    % White surface: world vertices + consistent vertex normals.
    sSurf = in_tess_bst(WhiteSurfFile);
    Vw = cs_convert(sMriRef,'scs','world', sSurf.Vertices);
    Nw = local_vertnormals(Vw, sSurf.Faces);

    % Tangential heat smoother (needs the FULL surface; no subsampling when active).
    HS = [];
    if Opts.HeatT > 0
        HS = local_heat_setup(WhiteSurfFile, Opts.HeatT);
    elseif size(Vw,1) > Opts.SubSample
        sel = round(linspace(1, size(Vw,1), Opts.SubSample));
        Vw = Vw(sel,:); Nw = Nw(sel,:);
    end
    d = Opts.SampleDist * wpm;

    S = struct();
    S.gmW = Vw + d*Nw; S.wmW = Vw - d*Nw;
    S.c = mean(Vw,1)'; S.wpm = wpm; S.Tinit = TinitWorld;
    S.sMriPet = sMriPet; S.PET = double(sMriPet.Cube); S.cubeSz = size(S.PET(:,:,:,1));
    S.HS = HS; S.FrameCombine = Opts.FrameCombine;

    % Contrast sign at the init, then maximize |contrast|.
    Q0 = local_Q([0 0 0 0 0 0], S);
    sgn = sign(mean(Q0(isfinite(Q0)))); if sgn==0, sgn=1; end
    costf = @(th) local_cost(local_Q(th, S), sgn, Opts.Mslope);
    costInit = costf([0 0 0 0 0 0]);

    [theta, cost] = local_optimize(costf);
    TransfWorld = local_T(theta, S.c, wpm, TinitWorld);
    info = struct('theta',theta,'costInit',costInit,'nVert',size(Vw,1),'sign',sgn,'cost',cost, ...
                  'heatT',Opts.HeatT,'nFrames',size(S.PET,4));
end


%% ===== optimizer: coordinate descent + expanding line search (Powell-like) =====
function [theta, cost] = local_optimize(costf)
    theta = zeros(1,6); cost = costf(theta);
    rotSteps = deg2rad([4 2 1 0.5]); trSteps = [8 4 2 1 0.5];
    for sweep = 1:4
        improved = false;
        for p = 1:6
            if p <= 3, steps = rotSteps; else, steps = trSteps; end
            for st = steps
                for dir = [1 -1]
                    while true
                        cand = theta; cand(p) = cand(p) + dir*st; cc = costf(cand);
                        if cc < cost - 1e-5, cost = cc; theta = cand; improved = true; else, break; end
                    end
                end
            end
        end
        if ~improved, break; end
    end
end


%% ===== cost / sampling =====
function Q = local_Q(theta, S)
    Tinv = inv(local_T(theta, S.c, S.wpm, S.Tinit));
    gVox = cs_convert(S.sMriPet,'world','voxel', local_apply(Tinv, S.gmW));
    wVox = cs_convert(S.sMriPet,'world','voxel', local_apply(Tinv, S.wmW));
    nF = size(S.PET,4); nP = size(S.gmW,1);
    Ig = nan(nP,nF); Iw = nan(nP,nF);
    for f = 1:nF
        Ig(:,f) = local_interp(S.PET(:,:,:,f), gVox, S.cubeSz);
        Iw(:,f) = local_interp(S.PET(:,:,:,f), wVox, S.cubeSz);
    end
    if strcmpi(S.FrameCombine,'median') && nF > 1
        if ~isempty(S.HS), Ig = local_heatsmooth(S.HS, Ig); Iw = local_heatsmooth(S.HS, Iw); end
        Qf = 100*(Ig - Iw) ./ (0.5*(Ig + Iw) + eps);
        Q = median(Qf, 2, 'omitnan');
    else
        ig = mean(Ig,2,'omitnan'); iw = mean(Iw,2,'omitnan');
        if ~isempty(S.HS), ig = local_heatsmooth(S.HS, ig); iw = local_heatsmooth(S.HS, iw); end
        Q = 100*(ig - iw) ./ (0.5*(ig + iw) + eps);
    end
end

function c = local_cost(Q, sgn, M)
    Q = Q(isfinite(Q));
    if numel(Q) < 100, c = 2; return; end
    c = mean(1 - tanh(sgn * M * Q));
end


%% ===== geometry helpers =====
function T = local_T(theta, c, wpm, Tinit)
    rx=theta(1); ry=theta(2); rz=theta(3); t=theta(4:6)' * wpm;
    Rx=[1 0 0;0 cos(rx) -sin(rx);0 sin(rx) cos(rx)];
    Ry=[cos(ry) 0 sin(ry);0 1 0;-sin(ry) 0 cos(ry)];
    Rz=[cos(rz) -sin(rz) 0;sin(rz) cos(rz) 0;0 0 1];
    Dl=[Rz*Ry*Rx, t; 0 0 0 1]; Tc=[eye(3) c; 0 0 0 1];
    T = Tc*Dl/Tc * Tinit;
end

function P2 = local_apply(M, P)
    P2 = (M * [P'; ones(1,size(P,1))])'; P2 = P2(:,1:3);
end

function I = local_interp(PET, vox, cubeSz)
    ok = all(vox>=1,2) & vox(:,1)<=cubeSz(1) & vox(:,2)<=cubeSz(2) & vox(:,3)<=cubeSz(3);
    I = nan(size(vox,1),1);
    if any(ok), I(ok) = interpn(PET, vox(ok,1), vox(ok,2), vox(ok,3), 'linear', NaN); end
end

function N = local_vertnormals(V, F)
    % cross(v3-v1, v2-v1) matches Brainstorm's face winding -> OUTWARD normals (verified
    % against the stored VertNormals and the white->pial direction). The opposite order
    % gives inward normals.
    v1=V(F(:,1),:); v2=V(F(:,2),:); v3=V(F(:,3),:);
    fn = cross(v3-v1, v2-v1, 2); nV = size(V,1); N = zeros(nV,3); idx = F(:);
    for dd=1:3, N(:,dd) = accumarray(idx, repmat(fn(:,dd),3,1), [nV 1]); end
    N = N ./ (sqrt(sum(N.^2,2))+eps);
end


%% ===== tangential (Laplace-Beltrami heat) smoothing of a per-vertex field =====
function HS = local_heat_setup(WhiteSurfFile, t)
% Build the cached implicit-Euler heat operator (M + t*K) per hemisphere from the white
% surface LBO. Smoothing diffuses a vertex field ALONG the surface only.
    LBO = tess_operators(WhiteSurfFile, 'Laplace-Beltrami');
    HS.gv = LBO.GlobalVertices; HS.M = LBO.Mass; HS.dA = cell(1,numel(HS.M));
    for hh = 1:numel(HS.M)
        HS.dA{hh} = decomposition(HS.M{hh} + t*LBO.Operator{hh}, 'chol');
    end
end

function g = local_heatsmooth(HS, f)
    f(~isfinite(f)) = 0;                          % heat solve needs finite RHS
    g = f;
    for hh = 1:numel(HS.M)
        vH = double(HS.gv{hh}(:));
        g(vH,:) = HS.dA{hh} \ (HS.M{hh} * f(vH,:));
    end
end
