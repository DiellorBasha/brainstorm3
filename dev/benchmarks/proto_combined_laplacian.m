function R = proto_combined_laplacian(trkFile, seedRegion, beta)
% PROTO_COMBINED_LAPLACIAN: build a structural connectome from a tractography .trk and form the
% COMBINED operator L = L_LBO + connectome, then demo "inject heat -> spreads locally (geodesic) AND
% to connectome-target regions (network)". The substrate is symmetric (LBO + connectome graph
% Laplacian), so directional spread emerges from the seed, not the operator.
%
% Pipeline: parse .trk endpoints -> voxel->MNI (DSI ICBM152 2mm affine) -> subject SCS -> nearest
% Desikan scout -> region x region connectome C (on the ICBM152 default subject). LBO + region->vertex
% lifting of C on a MANIFOLD cortex (sub-MTL0002). Split heat: implicit LBO + explicit connectome mix.
%
% USAGE: R = proto_combined_laplacian('.../whole_brain.trk', 'superiorfrontal L', 0.30)
%
% Author: Diellor Basha, 2026 (prototype; ICBM152 template connectome)
    if nargin<2||isempty(seedRegion), seedRegion='superiorfrontal L'; end
    if nargin<3||isempty(beta), beta=0.30; end
    here=bst_fileparts(mfilename('fullpath'));

    % --- connectome on the ICBM152 default subject (Desikan, aligned to the MNI tractography) ---
    [sDef,~]=bst_get('Subject',0); sMri=in_mri_bst(sDef.Anatomy(sDef.iAnatomy).FileName);   % default anatomy (holds the ICBM152 + HCP-1065 fibers)
    sIc=in_tess_bst(sDef.Surface(find(~cellfun('isempty',regexp({sDef.Surface.FileName},'cortex_mid_low\.mat$','once')),1)).FileName);
    EP=local_read_trk_endpoints(trkFile);
    toMNI=@(P)[-P(:,1)+79.5, -P(:,2)+81.5, P(:,3)-72];        % DSI ICBM152 2mm voxel(*2=mm)->MNI mm
    scs1=cs_convert(sMri,'mni','scs',toMNI(EP(:,1:3))/1000);
    scs2=cs_convert(sMri,'mni','scs',toMNI(EP(:,4:6))/1000);
    ai=find(strcmp({sIc.Atlas.Name},'Desikan-Killiany'),1); Sc=sIc.Atlas(ai).Scouts; nReg=numel(Sc);
    vReg=zeros(size(sIc.Vertices,1),1); for r=1:nReg, vReg(Sc(r).Vertices)=r; end
    r1=vReg(dsearchn(sIc.Vertices,scs1)); r2=vReg(dsearchn(sIc.Vertices,scs2));
    ok=r1>0 & r2>0 & r1~=r2; C=accumarray([r1(ok) r2(ok)],1,[nReg nReg]); C=C+C';

    % --- combined operator on the MANIFOLD MTL0002 cortex (LBO works; Desikan matched by label) ---
    [sM,~]=bst_get('Subject','sub-MTL0002');
    wf=sM.Surface(find(~cellfun('isempty',regexp({sM.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    sWm=in_tess_bst(wf); nV=size(sWm.Vertices,1);
    aiM=find(strcmp({sWm.Atlas.Name},'Desikan-Killiany'),1); ScM=sWm.Atlas(aiM).Scouts;
    P=sparse(nV,nReg);
    for r=1:nReg, j=find(strcmp({ScM.Label},Sc(r).Label),1); if ~isempty(j), P(ScM(j).Vertices,r)=1; end; end
    LBO=tess_operators(wf,'Laplace-Beltrami'); M=sparse(nV,nV); K=sparse(nV,nV);
    for hh=1:2, gv=double(LBO.GlobalVertices{hh}(:)); M(gv,gv)=LBO.Mass{hh}; K(gv,gv)=LBO.Operator{hh}; end
    Crw=C./max(sum(C,2),eps); Pn=P./max(sum(P,1),eps);

    % --- inject-heat demo: LBO-only vs combined ---
    sidx=find(strcmp({Sc.Label},seedRegion),1); u0=full(double(P(:,sidx)>0)); u0=u0/sum(u0);
    dt=4e-5; nStep=8; dLBO=decomposition(M+dt*K); uL=u0; uC=u0;
    for n=1:nStep
        uL=dLBO\(M*uL);
        uC=dLBO\(M*uC); uC=uC + beta*(P*(Crw*(Pn'*uC)) - uC);
    end

    f=figure('Visible','off','Position',[40 40 1000 460]); cmax=max(uC); V=sWm.Vertices; Fa=sWm.Faces;
    ttl={'LBO only (geodesic) - stays in one hemisphere','Combined LBO+connectome - jumps via the network'};
    for k=1:2
        subplot(1,2,k); u=uL*(k==1)+uC*(k==2);
        patch('Faces',Fa,'Vertices',V,'FaceVertexCData',u,'FaceColor','interp','EdgeColor','none');
        clim([0 cmax]); colormap(hot); view([0 90]); axis equal off vis3d; camlight headlight; lighting gouraud; title(ttl{k});
    end
    png=fullfile(here,'proto_combined_laplacian.png'); print(f,png,'-dpng','-r110'); close(f);
    fprintf('connectome %dx%d (density %.0f%%) | seed %s | figure -> %s\n', nReg, nReg, 100*nnz(triu(C,1))/(nReg*(nReg-1)/2), seedRegion, png);
    R=struct('C',C,'Sc',{Sc},'uL',uL,'uC',uC,'png',png);
end

function EP=local_read_trk_endpoints(trk)
    fid=fopen(trk,'r','l'); fseek(fid,36,'bof'); nSc=fread(fid,1,'int16');
    fseek(fid,238,'bof'); nPr=fread(fid,1,'int16'); fseek(fid,988,'bof'); nT=fread(fid,1,'int32');
    fseek(fid,1000,'bof'); EP=zeros(nT,6);
    for i=1:nT
        m=fread(fid,1,'int32'); pts=reshape(fread(fid,(3+nSc)*m,'float32'),3+nSc,m)';
        if nPr>0, fread(fid,nPr,'float32'); end
        EP(i,:)=[pts(1,1:3) pts(end,1:3)];
    end
    fclose(fid);
end
