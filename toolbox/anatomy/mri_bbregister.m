function [TransfWorld, cost, info] = mri_bbregister(sMriPet, sMriRef, WhiteSurfFile, TinitWorld, Opts)
% MRI_BBREGISTER: Boundary-Based Registration of a (PET) volume to anatomy (Greve 2009).
%
% Native replication of FreeSurfer bbregister: finds the 6-DOF RIGID transform that aligns
% the moving volume (PET) to the subject's anatomy by MAXIMIZING the GM/WM intensity contrast
% sampled across the white surface - instead of a global intensity similarity. This is robust
% for low-contrast PET (e.g. flortaucipir/tau) where intensity-based coregistration drifts.
%
% It is a LOCAL refiner: pass a coarse initial alignment (e.g. from SPM coregistration) and it
% polishes it. Multi-start guards against a poor init. Rotations/translations are optimized in
% radians/millimetres; all geometry is done in Brainstorm 'world' coordinates.
%
% USAGE:  [TransfWorld, cost, info] = mri_bbregister(sMriPet, sMriRef, WhiteSurfFile, TinitWorld, Opts)
%
% INPUTS:
%   sMriPet      : moving MRI/PET struct (Cube + SCS/world geometry).
%   sMriRef      : reference anatomical MRI struct (defines SCS<->world for the surface).
%   WhiteSurfFile: white (GM/WM) surface file, registered to sMriRef.
%   TinitWorld   : 4x4 initial PET-world -> anat-world transform ([] = identity, i.e. the PET
%                  is already approximately in anatomical space).
%   Opts         : (optional) .SampleDist (mm, default 1.5), .SubSample (max verts, 5000),
%                  .MultiStart (n, default 4), .Mslope (default 0.5).
%
% OUTPUTS:
%   TransfWorld  : 4x4 refined PET-world -> anat-world rigid transform.
%   cost         : final BBR cost (0 = sharp boundary, 1 = no contrast).
%   info         : struct with .theta .costInit .nVert .sign .cost.
%
% Reference: Greve DN, Fischl B. Accurate and robust brain image alignment using boundary-based
%            registration. NeuroImage 2009;48(1):63-72.
%
% SEE ALSO: mri_coregister, cs_convert
%
% Author: Diellor Basha, 2026

    if (nargin < 4), TinitWorld = []; end
    if isempty(TinitWorld), TinitWorld = eye(4); end
    if (nargin < 5) || isempty(Opts), Opts = struct(); end
    D = struct('SampleDist',1.5, 'SubSample',5000, 'MultiStart',4, 'Mslope',0.5);
    fn = fieldnames(D); for i=1:numel(fn), if ~isfield(Opts,fn{i})||isempty(Opts.(fn{i})), Opts.(fn{i})=D.(fn{i}); end; end

    % World units per millimetre (Brainstorm 'world' may be metres): length of 1 mm in world.
    p0 = cs_convert(sMriRef, 'voxel', 'world', [1 1 1]);
    p1 = cs_convert(sMriRef, 'voxel', 'world', [2 1 1]);
    wpm = norm(p1 - p0) / sMriRef.Voxsize(1);                % world units / mm

    % ----- white surface: world-space vertices + consistent vertex normals -----
    sSurf = in_tess_bst(WhiteSurfFile);
    Vw = cs_convert(sMriRef, 'scs', 'world', sSurf.Vertices);
    Nw = local_vertnormals(Vw, sSurf.Faces);                 % unit, consistent winding
    nV = size(Vw,1);
    if nV > Opts.SubSample
        sel = round(linspace(1, nV, Opts.SubSample));
        Vw = Vw(sel,:); Nw = Nw(sel,:);
    end
    d = Opts.SampleDist * wpm;                               % step in world units
    gmW = Vw + d*Nw;                                         % OUTSIDE the surface (GM side, up to sign)
    wmW = Vw - d*Nw;                                         % INSIDE (WM side)

    PET = double(sMriPet.Cube(:,:,:,1));
    cubeSz = size(PET);
    c = mean(Vw,1)';                                         % rotate about the surface centroid (world)

    % ----- cost (estimate contrast sign at the init, then maximize |contrast|) -----
    Qfun = @(th) local_Q(th, gmW, wmW, c, wpm, TinitWorld, sMriPet, PET, cubeSz);
    Q0 = Qfun([0 0 0 0 0 0]);
    sgn = sign(mean(Q0(isfinite(Q0)))); if sgn==0, sgn=1; end
    costf = @(th) local_cost(Qfun(th), sgn, Opts.Mslope);
    costInit = costf([0 0 0 0 0 0]);

    % ----- optimize: coordinate descent with expanding line search (Powell-like, as in
    %       Greve's BBR). Robust where fminsearch's zero-init simplex stalls, and the
    %       expansion gives a wide capture range (recovers large mis-inits). -----
    [theta, cost] = local_optimize(costf);

    TransfWorld = local_T(theta, c, wpm, TinitWorld);
    info = struct('theta',theta, 'costInit',costInit, 'nVert',size(Vw,1), 'sign',sgn, 'cost',cost);
end


function [theta, cost] = local_optimize(costf)
% Coordinate descent over [rx ry rz (rad)  tx ty tz (mm)] with multi-scale expanding line
% search per parameter. Steps in mm for translations, degrees for rotations.
    theta = zeros(1,6); cost = costf(theta);
    rotSteps = deg2rad([4 2 1 0.5]); trSteps = [8 4 2 1 0.5];
    for sweep = 1:4
        improved = false;
        for p = 1:6
            if p <= 3, steps = rotSteps; else, steps = trSteps; end
            for st = steps
                for dir = [1 -1]
                    while true
                        cand = theta; cand(p) = cand(p) + dir*st;
                        cc = costf(cand);
                        if cc < cost - 1e-5
                            cost = cc; theta = cand; improved = true;
                        else
                            break;
                        end
                    end
                end
            end
        end
        if ~improved, break; end
    end
end


%% ===== helpers =====
function N = local_vertnormals(V, F)
    v1=V(F(:,1),:); v2=V(F(:,2),:); v3=V(F(:,3),:);
    fn = cross(v2-v1, v3-v1, 2);                              % face normals (consistent winding)
    nV = size(V,1); N = zeros(nV,3); idx = F(:);
    for dd=1:3, N(:,dd) = accumarray(idx, repmat(fn(:,dd),3,1), [nV 1]); end
    N = N ./ (sqrt(sum(N.^2,2))+eps);
end

function T = local_T(theta, c, wpm, Tinit)
    rx=theta(1); ry=theta(2); rz=theta(3); t=theta(4:6)' * wpm;   % translation mm -> world
    Rx=[1 0 0;0 cos(rx) -sin(rx);0 sin(rx) cos(rx)];
    Ry=[cos(ry) 0 sin(ry);0 1 0;-sin(ry) 0 cos(ry)];
    Rz=[cos(rz) -sin(rz) 0;sin(rz) cos(rz) 0;0 0 1];
    Dl=[Rz*Ry*Rx, t; 0 0 0 1];
    Tc=[eye(3) c; 0 0 0 1];
    T = Tc*Dl/Tc * Tinit;                                    % perturb about centroid, then init
end

function Q = local_Q(theta, gmW, wmW, c, wpm, Tinit, sMriPet, PET, cubeSz)
    Tinv = inv(local_T(theta, c, wpm, Tinit));               % anat-world -> pet-world
    gV = cs_convert(sMriPet,'world','voxel', local_apply(Tinv, gmW));   % 1-based voxel indices
    wV = cs_convert(sMriPet,'world','voxel', local_apply(Tinv, wmW));
    Ig = local_interp(PET, gV, cubeSz);
    Iw = local_interp(PET, wV, cubeSz);
    Q = 100*(Ig - Iw) ./ (0.5*(Ig + Iw) + eps);
end

function P2 = local_apply(M, P)
    P2 = (M * [P'; ones(1,size(P,1))])'; P2 = P2(:,1:3);
end

function I = local_interp(PET, vox, cubeSz)
    ok = all(vox>=1,2) & vox(:,1)<=cubeSz(1) & vox(:,2)<=cubeSz(2) & vox(:,3)<=cubeSz(3);
    I = nan(size(vox,1),1);
    if any(ok)
        I(ok) = interpn(PET, vox(ok,1), vox(ok,2), vox(ok,3), 'linear', NaN);
    end
end

function c = local_cost(Q, sgn, M)
    Q = Q(isfinite(Q));
    if numel(Q) < 100, c = 2; return; end
    c = mean(1 - tanh(sgn * M * Q));                         % minimized when |contrast| large
end
