function [GMpvc, info] = pet_pvc_surface(sMriPet, sMriRef, WhiteSurfFile, Opts)
% PET_PVC_SURFACE: surface-native partial volume correction in the intrinsic/extrinsic split.
%
% CANDIDATE / experimental. Corrects partial-volume effect directly on the cortical surface.
% The TISSUE-FRACTION effect (spill-in/out) is, by construction, a purely RADIAL/extrinsic
% phenomenon: along the cortical normal the PSF crosses WM->GM->CSF, while the tangential plane
% is all the same tissue (its PSF integral is unity -> zero spill). So the partial-volume
% correction is the per-vertex 1-D Mueller-Gaertner:
%     GM = (obs - m_wm*WM)/m_gm ,   m_wm,m_gm = PSF tissue fractions from THICKNESS (|pial-white|).
% This is theoretically exact for locally-flat cortex; curvature couples the directions only where
% the radius of curvature shrinks to ~a PSF-width of the thickness (a small high-curvature minority).
%
% Two things the experiments established and that this first-order version bakes in:
%   - MEDIAL WALL: vertices where FreeSurfer collapses pial onto white (thickness ~ 0) are NOT
%     cortex (verified to coincide 100% with the Desikan-unlabeled region). They are MASKED (NaN),
%     not regularized.
%   - TANGENTIAL deconvolution is a SEPARATE concern (sharpening real signal, not tissue-fraction
%     correction) and is OFF by default. It is also ill-posed and scale-mismatched to global
%     eigenmodes; enable Opts.DoTangential only for optional sub-PSF sharpening.
%
% Sampling is anchored to the white->pial axis (NOT vertex normals), so tissue identity is robust
% to surface conventions. Radial step needs no eigenmodes -> scales to the high-res surface.
%
% USAGE:  [GMpvc, info] = pet_pvc_surface(sMriPet, sMriRef, WhiteSurfFile, Opts)
%
% INPUTS:
%   sMriPet      : registered PET MRI struct (Cube + geometry; .PET metadata used for PSF if present).
%   sMriRef      : anatomical MRI struct (defines surface SCS<->world).
%   WhiteSurfFile: white surface file (pial derived by name unless Opts.PialSurfFile given).
%   Opts         : .PialSurfFile ('') .PsfFwhm (mm, [] = auto from scanner/metadata)
%                  .WmOffset (mm into WM for the WM-level estimate, 3) .MinThick (mm, medial-wall
%                  mask below this, 1.0) .GmFracMin (radial-division guard, 0.1)
%                  .DoTangential (false) .Alpha (Wiener reg, 1e-2) .nModes (600).
%
% OUTPUTS:
%   GMpvc        : [nVert x 1] partial-volume-corrected GM uptake on the cortex (pre-SUVR);
%                  NaN at masked medial-wall / non-cortex vertices.
%   info         : struct(.mask, .nMask, .thin, .WmLevel, .PsfFwhm, .m_gm, .m_wm, .thickness, .observed).
%
% SEE ALSO: pet_pvc, pet_gtm, tess_operators, tess_eigen, mri_bbregister
%
% Author: Diellor Basha, 2026

    if (nargin<4)||isempty(Opts), Opts=struct(); end
    Def=struct('PialSurfFile','','PsfFwhm',[],'WmOffset',3,'MinThick',1.0,'GmFracMin',0.1, ...
               'DoTangential',false,'Alpha',1e-2,'nModes',600);
    fn=fieldnames(Def); for i=1:numel(fn), if ~isfield(Opts,fn{i})||isempty(Opts.(fn{i})), Opts.(fn{i})=Def.(fn{i}); end; end

    % World units per mm.
    wpm = norm(cs_convert(sMriRef,'voxel','world',[2 1 1])-cs_convert(sMriRef,'voxel','world',[1 1 1]))/sMriRef.Voxsize(1);

    % Geometry: white, pial -> mid; white->pial axis (GM-ward) + thickness.
    sW=in_tess_bst(WhiteSurfFile); Vw=cs_convert(sMriRef,'scs','world', sW.Vertices);
    PialSurfFile=Opts.PialSurfFile; if isempty(PialSurfFile), PialSurfFile=strrep(WhiteSurfFile,'white','pial'); end
    sP=in_tess_bst(PialSurfFile); Vp=cs_convert(sMriRef,'scs','world', sP.Vertices);
    Uvec=Vp-Vw; th=sqrt(sum(Uvec.^2,2)); u=Uvec./max(th,eps); Xmid=0.5*(Vw+Vp);

    % PSF width (auto from scanner/metadata if not given).
    fwhm=Opts.PsfFwhm;
    if isempty(fwhm)
        if isfield(sMriPet,'PET'), fwhm=pet_scanner_fwhm(sMriPet.PET, 6); else, fwhm=6; end
    end
    sig=(fwhm/2.355)*wpm;                                   % PSF sigma in world units

    PET=double(sMriPet.Cube(:,:,:,1)); cs=size(PET);
    smp=@(P) local_interp(PET, cs_convert(sMriPet,'world','voxel',P), cs);

    % Observed PET at mid surface; WM level from a robust sample well inside white.
    observed = smp(Xmid);
    wmSamp = smp(Vw - Opts.WmOffset*wpm.*u);
    WmLevel = median(wmSamp(isfinite(wmSamp) & wmSamp>0));

    % Radial PSF tissue fractions from thickness (1-D Gaussian along the normal; CSF symmetric).
    m_wm = 0.5*erfc((th/2)./(sig*sqrt(2)));
    m_gm = 1 - 2*m_wm;
    thin = m_gm < Opts.GmFracMin;                          % residual ill-conditioning (rare after mask)
    m_gm_reg = max(m_gm, Opts.GmFracMin);

    % Radial Mueller-Gaertner: remove WM spill-in, recover GM spill-out.
    GMpvc = (observed - m_wm.*WmLevel) ./ m_gm_reg;

    % Medial wall / non-cortex (FreeSurfer collapses pial onto white -> thickness ~ 0; coincides
    % with the Desikan-unlabeled region). Not cortex -> masked to NaN below (after any tangential step).
    mask = (th/wpm) < Opts.MinThick;

    % OPTIONAL tangential intrinsic deconvolution (regularized LBO Wiener), per hemisphere.
    % Separate from the tissue-fraction PVC; off by default.
    if Opts.DoTangential
        Eig=tess_eigen(WhiteSurfFile,'Laplace-Beltrami','nModes',Opts.nModes);
        LBO=tess_operators(WhiteSurfFile,'Laplace-Beltrami');
        tTan=sig^2/2;
        for hh=1:numel(Eig.Phi)
            gv=double(Eig.GlobalVertices{hh}(:)); Phi=Eig.Phi{hh}; lam=Eig.Lambda{hh}(:); Mh=LBO.Mass{hh};
            G=exp(-tTan*lam); h=G./(G.^2+Opts.Alpha);
            f=GMpvc(gv); f(~isfinite(f))=0;
            GMpvc(gv)=Phi*(h.*(Phi'*(Mh*f)));
        end
    end
    GMpvc(mask) = NaN;

    info=struct('mask',mask,'nMask',nnz(mask),'thin',thin & ~mask,'WmLevel',WmLevel,'PsfFwhm',fwhm, ...
                'm_gm',m_gm,'m_wm',m_wm,'thickness',th/wpm,'observed',observed,'nThin',nnz(thin & ~mask));
end

function I=local_interp(PET,vox,cs)
    ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3); I=nan(size(vox,1),1);
    if any(ok), I(ok)=interpn(PET,vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN); end
end
