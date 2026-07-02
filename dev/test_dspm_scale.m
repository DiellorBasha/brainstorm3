function tests = test_dspm_scale
tests = functiontests(localfunctions);
end
function test_anchor_direction_and_scale(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig,'DynamicsOverlay');
    src = panel_bst_dynamics('i_src_resultfile', D);
    R = in_bst_results(src,0,'ImagingKernel','ImagingKernelMode','Eigenvalues','ModeHemisphere','DiracEigenFile');
    E = in_bst_eigen(R.DiracEigenFile);
    lam=double(R.Eigenvalues(:)); hemi=double(R.ModeHemisphere(:));
    nV=0; for h=1:2, nV=max(nV,max(E.GlobalVertices{h}(:))); end
    Km=double(R.ImagingKernelMode); Krec=zeros(3*nV,size(Km,2));
    for h=1:2
        ord=find(hemi==h);[~,s]=sort(lam(ord),'ascend');ord=ord(s);
        Ph=double(E.Phi{h}); gv=E.GlobalVertices{h}(:); Uf=Ph*Km(ord,:);
        Krec((gv-1)*3+1,:)=Uf(2:4:end,:);Krec((gv-1)*3+2,:)=Uf(3:4:end,:);Krec((gv-1)*3+3,:)=Uf(4:4:end,:);
    end
    Kv=double(R.ImagingKernel);
    % direction identical per vertex
    cosv=zeros(nV,1);
    for v=1:nV, a=Krec((v-1)*3+(1:3),:); b=Kv((v-1)*3+(1:3),:); na=norm(a,'fro'); nb=norm(b,'fro');
        if na>0&&nb>0, cosv(v)=sum(a(:).*b(:))/(na*nb); end, end
    verifyGreaterThan(t, median(cosv(cosv~=0)), 1-1e-6);
    % i_dspm_scale reproduces the per-vertex ratio: SIR .* recon == ImagingKernel
    sir = panel_bst_dynamics('i_dspm_scale', st, D);
    Kfix = Krec .* reshape(repmat(sir(:)',3,1),[],1);
    verifyLessThan(t, norm(Kfix(:)-Kv(:))/norm(Kv(:)), 1e-6);
end
