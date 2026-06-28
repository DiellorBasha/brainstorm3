function Rp = pseudo_longitudinal_check(subjects)
% PSEUDO_LONGITUDINAL_CHECK: order subjects by global cortical SUVR severity to form a PSEUDO-TIME
% series of LH Abeta & tau surface maps (sampled per subject, mapped to a manifold REFERENCE subject
% via registration spheres), then run the reaction-diffusion coupling inverter as a real-data
% sanity cross-check.
%
% IMPORTANT: this is cross-sectional pseudo-time, with few non-uniformly-spaced points - i.e. exactly
% the sparse/coarse regime that breaks the dense-time inverter (see validate_pet_spread_coupling).
% The recovered kappa is a pipeline demonstration on real maps, NOT a trustworthy estimate.
%
% USAGE: Rp = pseudo_longitudinal_check()                 % whole cohort
%        Rp = pseudo_longitudinal_check({'sub-MTL0166',...})
%
% Author: Diellor Basha, 2026
    ref='sub-MTL0002';                                  % manifold reference (the template is non-manifold)
    [sR,~]=bst_get('Subject',ref);
    wfR=sR.Surface(find(~cellfun('isempty',regexp({sR.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    LBO=tess_operators(wfR,'Laplace-Beltrami'); gv=double(LBO.GlobalVertices{1}(:));
    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name};
        subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once')));
    end
    A=[]; T=[]; sev=[];
    for s=1:numel(subjects)
        m=local_surf_suvr_ref(subjects{s},'PET 18FNAV4694_mean_pvc', wfR, gv);
        n=local_surf_suvr_ref(subjects{s},'PET 18Fflortaucipir_mean_pvc', wfR, gv);
        if isempty(m)||isempty(n), continue; end
        A=[A,m]; T=[T,n]; sev=[sev,mean(m,'omitnan')]; %#ok<AGROW>
    end
    [~,ord]=sort(sev); A=A(:,ord); T=T(:,ord);          % pseudo-time = increasing amyloid severity
    if max(A(:))>0, A=A/max(A(:)); end; if max(T(:))>0, T=T/max(T(:)); end
    A=min(max(A,0),1); T=min(max(T,0),1);
    est=pet_spread_invert(wfR, A, T, 1);
    fprintf('\npseudo-longitudinal (%d subjects, severity-ordered, PSEUDO-TIME): kappa=%.2f rt=%.3f Dt=%.2e\n', size(A,2), est.kappa, est.rt, est.Dt);
    fprintf('  CAVEAT: cross-sectional pseudo-time + few non-uniform points = the sparse regime that\n');
    fprintf('          breaks this dense-time inverter; treat kappa as a pipeline demo, not an estimate.\n');
    Rp=struct('kappa',est.kappa,'est',est,'nSub',size(A,2),'A',A,'T',T);
end

function m=local_surf_suvr_ref(subj, petC, wfR, gv)
    m=[]; [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment};
    if ~any(strcmp(cmt,petC))||~any(strcmp(cmt,'ASEG')), return; end
    af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;
    sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
    wfS=sf('cortex_white_low\.mat$');
    sW=in_tess_bst(wfS); sP=in_tess_bst(sf('cortex_pial_low\.mat$'));
    sSuvr=pet_suvr(in_mri_bst(af(petC)), in_mri_bst(af('ASEG')));
    Vmid=0.5*(sW.Vertices+sP.Vertices); vox=cs_convert(sSuvr,'scs','voxel',Vmid);
    cs=size(sSuvr.Cube); ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3);
    v=nan(size(Vmid,1),1); v(ok)=interpn(double(sSuvr.Cube),vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN); v(~isfinite(v))=0;
    Wm=tess_interp_tess2tess(wfS, wfR, 0, 0, 0);         % map subject -> reference surface (reg spheres)
    vR=Wm*v; m=vR(gv); m(~isfinite(m))=0;
end
