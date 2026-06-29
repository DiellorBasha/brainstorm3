function FibFile = build_hcp1065_default_fibers(trkFile, nPts)
% BUILD_HCP1065_DEFAULT_FIBERS: bundle the HCP-1065 template tractography as a fibers asset on the
% protocol's DEFAULT anatomy (@default_subject), so it ships/loads like the average cortex and atlases.
% Fibers are anatomy (Brainstorm tracks them in the subject Fibers field) - this makes the structural
% connectome substrate available to derive connectome/combined Laplacians on demand (tess_connectome),
% and is what lets spectral/differential operators become truly WHOLE-BRAIN (the per-hemisphere LBO +
% the interhemispheric fibers).
%
% Alignment: DSI Studio exports the .trk in ICBM152 2mm voxel space; the validated transform is
%   MNI_mm = [-x+79.5, -y+81.5, z-72]  (voxel*2=mm folded in), then /1000 -> MNI metres -> cs_convert
%   'mni'->'scs'. Verified: fiber SCS X[-81 99] Y[-68 68] overlaps default cortex X[-80 99] Y[-71 69].
% (Brainstorm's import_fibers CS='trk'/'mni' path misaligned this export, hence the direct build.)
%
% Registration: db_add_surface uses iSubject==0 for the default subject (NOT the resolved index of a
% subject that merely uses default anatomy - passing that index triggers a db_surface_default index
% error). Manual fallback appends to ProtocolSubjects.DefaultSubject + sets iFibers.
%
% USAGE: build_hcp1065_default_fibers('/.../hcp1065.trk', 40)
%
% Author: Diellor Basha, 2026
    if nargin<2 || isempty(nPts), nPts=40; end
    sMri = in_mri_bst(bst_get('Subject',0).Anatomy(bst_get('Subject',0).iAnatomy).FileName);

    % --- read full streamlines (raw .trk coords), resample to nPts ---
    fid=fopen(trkFile,'r','l');
    fseek(fid,6,'bof');  dimv=fread(fid,3,'int16')'; fseek(fid,12,'bof'); vox=fread(fid,3,'float32')';
    fseek(fid,36,'bof'); nSc=fread(fid,1,'int16');   fseek(fid,238,'bof'); nPr=fread(fid,1,'int16');
    fseek(fid,988,'bof'); nT=fread(fid,1,'int32');   fseek(fid,1000,'bof');
    Pts=zeros(nT,nPts,3,'single');
    for i=1:nT
        m=fread(fid,1,'int32'); raw=reshape(fread(fid,(3+nSc)*m,'float32'),3+nSc,m)';
        if nPr>0, fread(fid,nPr,'float32'); end
        p=raw(:,1:3); if m<2, p=[p;p]; m=size(p,1); end
        Pts(i,:,:)=interp1(linspace(0,1,m),p,linspace(0,1,nPts));
    end
    fclose(fid);

    % --- validated transform to SCS ---
    P2=double(reshape(Pts,[],3)); toMNI=@(P)[-P(:,1)+79.5,-P(:,2)+81.5,P(:,3)-72];
    PtsS=single(reshape(cs_convert(sMri,'mni','scs',toMNI(P2)/1000),nT,nPts,3));

    % --- build + save fibers file in the default-subject anatomy folder ---
    FibMat=db_template('fibersmat');
    FibMat.Points=PtsS;
    FibMat.Header=struct('id_string','TRACK','dim',dimv,'voxel_size',vox,'n_scalars',nSc,'n_properties',nPr,'n_count',nT);
    FibMat.Comment='HCP-1065';
    FibMat=fibers_helper('ComputeColor',FibMat);
    FibMat=bst_history('add',FibMat,'import','HCP-1065 template tractography (DSI Studio, 1065 HCP subjects, ICBM152 2mm) - validated SCS alignment');
    pinfo=bst_get('ProtocolInfo');
    FibFile=file_unique(bst_fullfile(pinfo.SUBJECTS,'@default_subject','tess_fibers_hcp1065.mat'));
    bst_save(FibFile,FibMat,'v7');

    % --- register on the DEFAULT subject (iSubject==0); manual fallback ---
    relFib=file_short(FibFile); ok=false;
    try; db_add_surface(0,relFib,'HCP-1065'); ok=true; catch; end
    if ~ok
        PS=bst_get('ProtocolSubjects'); ns=db_template('Surface');
        ns.FileName=relFib; ns.Comment='HCP-1065'; ns.SurfaceType='Fibers';
        PS.DefaultSubject.Surface(end+1)=ns; PS.DefaultSubject.iFibers=numel(PS.DefaultSubject.Surface);
        bst_set('ProtocolSubjects',PS);
    end
    db_save();
    fprintf('HCP-1065 fibers (%d streamlines x %d pts) registered on default anatomy: %s\n', nT, nPts, relFib);
end
