function test_helmholtz_covariant()
% TEST_HELMHOLTZ_COVARIANT: bst_helmholtz decomposition on the 'Covariant' node.
%   Gate 2 (Div matches validated reference), Gate 4 (sign: source -> Div>0),
%   Gate 5 (residual div-free & curl-free = feature-complete), + output contract.
%
% USAGE:  test_helmholtz_covariant   % Brainstorm running, TutorialAuditory loaded
% Authors: Diellor Basha, 2026

    global GlobalData;
    PF={'FAIL','PASS'}; pass=true;
    relData='Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat'; sStudy=bst_get('DataFile',relData);
    R=''; for j=1:numel(sStudy.Result)
        if ~isempty(regexp(sStudy.Result(j).Comment,'MN: MEG\(Unconstr\)','once')) && ~isempty(regexp(sStudy.Result(j).FileName,'KERNEL','once'))
            R=sStudy.Result(j).FileName; break; end
    end
    if isempty(R), fprintf('SKIPPED (no unconstrained kernel)\n'); fprintf('\n==== SUITE: %s ====\n',PF{pass+1}); return; end
    [iDS,iRes]=bst_memory('LoadResultsFileFull',['link|' R '|' relData]);
    SurfaceFile=GlobalData.DataSet(iDS).Results(iRes).SurfaceFile;
    Cov=bst_get_operator_node(SurfaceFile,'Covariant'); LBO=bst_get_operator_node(SurfaceFile,'Laplace-Beltrami');
    Surf=in_tess_bst(SurfaceFile,0); Mani=tess_manifold(SurfaceFile);
    Op=bst_helmholtz('Prepare',{Cov,LBO},Mani,Surf,'Domain','vertex');
    [tv,~]=bst_memory('GetTimeVector',iDS,iRes,[]); [~,it]=min(abs(tv-202));
    Jall=double(bst_memory('GetResultsValues',iDS,iRes,[],it,0));
    Ht=bst_helmholtz('Frame',Op,Jall,false);

    % ---- GATE 2: Ht.Div matches the validated FEM divergence ----
    divFEM=zeros(size(Surf.Vertices,1),1);
    for hh=1:numel(Op.vH), vH=Op.vH{hh};
        divFEM(vH)=Op.Wfv{hh}*(Op.Gx{hh}*Jall(3*(vH-1)+1)+Op.Gy{hh}*Jall(3*(vH-1)+2)+Op.Gz{hh}*Jall(3*(vH-1)+3));
    end
    c2=corr(Ht.Div,divFEM); g2=c2>0.999;
    fprintf('GATE2 corr(Div, validated divFEM) = %.4f => %s\n', c2, PF{g2+1}); pass=pass&&g2;

    % ---- GATE 4: sign -- a radial-outward source gives Div > 0 ----
    Jsrc=zeros(3*size(Surf.Vertices,1),1); P=Surf.Vertices;
    for hh=1:numel(Op.vH), vH=Op.vH{hh}; seed=vH(round(numel(vH)/2));
        d=P(vH,:)-P(seed,:);                                  % radial outward (a source)
        Jsrc(3*(vH-1)+1)=d(:,1); Jsrc(3*(vH-1)+2)=d(:,2); Jsrc(3*(vH-1)+3)=d(:,3);
    end
    Hs=bst_helmholtz('Frame',Op,Jsrc,false);
    g4 = mean(Hs.Div(Hs.Div~=0)) > 0;
    fprintf('GATE4 sign: radial source mean(Div)=%+.3e => %s\n', mean(Hs.Div(Hs.Div~=0)), PF{g4+1}); pass=pass&&g4;

    % ---- GATE 5: residual div-free AND curl-free (feature-complete) ----
    gR=0; sR=0; nrm=0;
    for hh=1:numel(Op.vH), vH=Op.vH{hh}; W=spdiags(Op.Af{hh},0,numel(Op.Af{hh}),numel(Op.Af{hh}));
        Rf=[Op.Fvf{hh}*Ht.Hresid(vH,1), Op.Fvf{hh}*Ht.Hresid(vH,2), Op.Fvf{hh}*Ht.Hresid(vH,3)];
        gR=gR+norm(Op.Gx{hh}'*W*Rf(:,1)+Op.Gy{hh}'*W*Rf(:,2)+Op.Gz{hh}'*W*Rf(:,3))^2;
        sR=sR+norm(Op.Sx{hh}'*W*Rf(:,1)+Op.Sy{hh}'*W*Rf(:,2)+Op.Sz{hh}'*W*Rf(:,3))^2;
        nrm=nrm+norm(Ht.Vtot(vH,:),'fro')^2;
    end
    wd=sqrt(gR)/sqrt(nrm); wv=sqrt(sR)/sqrt(nrm); g5=(wd<1e-2)&&(wv<1e-2);
    fprintf('GATE5 residual weak-div=%.2e weak-vort=%.2e (feature-complete) => %s\n', wd, wv, PF{g5+1}); pass=pass&&g5;

    % ---- output contract ----
    g0 = isfield(Ht,'Jnormal') && isfield(Ht,'Hresid') && isfield(Ht,'HarmFrac') && ~isfield(Ht,'Vharm');
    hf = Ht.HarmFrac;  ok_hf = (hf>1e-3) && (hf<0.5);    % TRUE non-conforming residual, not the old ~0 lie
    fprintf('CONTRACT Jnormal+Hresid+HarmFrac present, Vharm gone = %d ; HarmFrac=%.3f (true, not ~0) = %d\n', g0, hf, ok_hf);
    pass = pass && g0 && ok_hf;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
