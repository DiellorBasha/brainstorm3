function test_kernel_comparison
% TEST_KERNEL_COMPARISON: Compare standard (depth-weighted) MNE, no-depth MNE, and
% Eigen-MNE imaging kernels in TutorialAuditory. Eigen-MNE uses a spectral prior
% instead of depth weighting, so the FAIR comparison is no-depth MNE vs Eigen-MNE
% (depth-weighted MNE kept for reference). Produces figures + a summary under
% dev/tests/results/ and asserts basic invariants.
%
% Conforms to the dev/tests pattern (function, addpath, brainstorm, assert,
% 'ALL TESTS PASSED'). Toolbox-free (uses corrcoef, 'omitnan', randperm).
% Adapted from the kernel-comparison template.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
outDir = fullfile(thisDir, 'results'); if ~exist(outDir,'dir'); mkdir(outDir); end

% ----- locate protocol + study holding both stored kernels -----
iProt = bst_get('Protocol','TutorialAuditory');
if isempty(iProt); disp('SKIP: TutorialAuditory not present.'); return; end
gui_brainstorm('SetCurrentProtocol', iProt);
sStudies = bst_get('ProtocolStudies'); T=[];
for iS=1:numel(sStudies.Study)
    s=sStudies.Study(iS); if isempty(s.Result)||isempty(s.HeadModel)||isempty(s.NoiseCov); continue; end
    rc=string({s.Result.Comment});
    if any(rc=="MN: MEG(Constr) 2018") && any(rc=="Eigen-MNE: MEG 2018"); T=iS; break; end
end
if isempty(T); disp('SKIP: study with both MNE and Eigen-MNE kernels not found.'); return; end
s=sStudies.Study(T); rc=string({s.Result.Comment});
iMN=find(rc=="MN: MEG(Constr) 2018"); iMN=iMN(end);
iEG=find(rc=="Eigen-MNE: MEG 2018"); iEG=iEG(end);
R_std=in_bst_results(s.Result(iMN).FileName,0);
R_eig=in_bst_results(s.Result(iEG).FileName,0);
K_std=R_std.ImagingKernel; K_eig=R_eig.ImagingKernel; iGood=R_std.GoodChannel;
assert(isequal(size(K_std),size(K_eig)),'kernel dims differ');
assert(R_std.nComponents==1,'expects constrained (nComponents=1)');
[nV,nCh]=size(K_std);

% ----- constrained good-channel leadfield + no-depth MNE -----
baseHmFile=R_std.HeadModelFile; HM=in_bst_headmodel(baseHmFile,0);
ncFile=s.NoiseCov(1).FileName; chFile=bst_get('ChannelFileForStudy',s.FileName);
NC=load(file_fullpath(ncFile)); ChannelMat=in_bst_channel(chFile);
HMg=HM; HMg.Gain=HM.Gain(iGood,:);
G=bst_gain_orient(double(HMg.Gain),HM.GridOrient,getfield_d(HM,'GridAtlas',[]));   % [nCh x nV]
OPT=bst_inverse_linear_2018();
OPT.InverseMethod='minnorm'; OPT.InverseMeasure='amplitude'; OPT.SourceOrient={'fixed'};
OPT.NoiseCovMat=struct('NoiseCov',NC.NoiseCov(iGood,iGood),'nSamples',[],'FourthMoment',[]);
OPT.ChannelTypes={ChannelMat.Channel(iGood).Type};
OPT.SnrMethod='fixed'; OPT.SnrFixed=3; OPT.UseDepth=0;
Rnd=bst_inverse_linear_2018(HMg,OPT); K_nd=Rnd.ImagingKernel;
assert(isequal(size(K_nd),size(K_std)) && all(isfinite(K_nd(:))),'no-depth kernel bad');
V=in_tess_bst(R_std.SurfaceFile); Vert=V.Vertices; Face=V.Faces;

col_std=[0.20 0.40 0.80]; col_nd=[0.10 0.55 0.30]; col_eig=[0.85 0.33 0.10];

% ================= FIG 1: SVD spectrum (log) =================
sv_std=svd(K_std); sv_nd=svd(K_nd); sv_eig=svd(K_eig);
en=@(s) cumsum(s.^2)/sum(s.^2);
r99=@(s) find(en(s)>=0.99,1,'first');
f1=figure('Color','w','Position',[100 100 1200 380]);
subplot(1,3,1); semilogy(sv_std,'Color',col_std,'LineWidth',1.5); hold on;
semilogy(sv_nd,'Color',col_nd,'LineWidth',1.5); semilogy(sv_eig,'Color',col_eig,'LineWidth',1.5);
grid on; xlabel('singular value index'); ylabel('singular value (log)');
legend({'MNE (depth)','MNE (no depth)','Eigen-MNE'},'Location','northeast'); title('SVD spectrum');
subplot(1,3,2); plot(en(sv_std),'Color',col_std,'LineWidth',1.5); hold on;
plot(en(sv_nd),'Color',col_nd,'LineWidth',1.5); plot(en(sv_eig),'Color',col_eig,'LineWidth',1.5);
yline(0.99,'k--'); grid on; xlabel('index'); ylabel('cumulative energy');
title(sprintf('eff. rank_{99}: %d / %d / %d', r99(sv_std), r99(sv_nd), r99(sv_eig)));
subplot(1,3,3); n=min([numel(sv_nd) numel(sv_eig)]);
semilogy(sv_eig(1:n)./sv_nd(1:n),'k','LineWidth',1.2); yline(1,'k--'); grid on;
xlabel('index'); ylabel('\sigma_{eig}/\sigma_{noDepth}'); title('SV ratio');
sgtitle('Imaging-kernel singular value spectra');
saveFig(f1,fullfile(outDir,'kernel_svd_spectrum'));

% ================= FIG 2: global summary (no-depth MNE vs Eigen-MNE) =================
f2=figure('Color','w','Position',[100 100 1200 400]);
subplot(1,3,1); x=K_nd(:); y=K_eig(:); m=numel(x);
idx=1:m; if m>60000; idx=randperm(m,60000); end
scatter(x(idx),y(idx),1,[.6 .6 .6],'filled','MarkerFaceAlpha',0.15); hold on;
ll=[min([x;y]) max([x;y])]; plot(ll,ll,'k--');
cc=corrcoef(x,y); grid on; xlabel('MNE (no depth) weights'); ylabel('Eigen-MNE weights');
title(sprintf('weight correlation r=%.3f', cc(1,2)));
subplot(1,3,2); fr=[norm(K_std,'fro') norm(K_nd,'fro') norm(K_eig,'fro')];
b=bar(fr,0.6); b.FaceColor='flat'; b.CData=[col_std;col_nd;col_eig];
set(gca,'XTickLabel',{'MNE depth','MNE noDepth','Eigen'}); ylabel('Frobenius norm'); title('global magnitude');
subplot(1,3,3); ed=linspace(ll(1),ll(2),200); ctr=0.5*(ed(1:end-1)+ed(2:end));
semilogy(ctr,histcounts(K_nd(:),ed,'Normalization','probability'),'Color',col_nd,'LineWidth',1.2); hold on;
semilogy(ctr,histcounts(K_eig(:),ed,'Normalization','probability'),'Color',col_eig,'LineWidth',1.2);
grid on; xlabel('kernel weight'); ylabel('probability (log)'); legend({'MNE noDepth','Eigen'}); title('weight distribution');
sgtitle('Global kernel comparison (no-depth MNE vs Eigen-MNE)');
saveFig(f2,fullfile(outDir,'kernel_global_summary'));

% ================= FIG 3: per-vertex agreement on cortex (#2) =================
cosab=zeros(nV,1); nr=zeros(nV,1);
nn1=sqrt(sum(K_nd.^2,2)); nn2=sqrt(sum(K_eig.^2,2));
dotab=sum(K_nd.*K_eig,2);
ok=nn1>0 & nn2>0; cosab(ok)=dotab(ok)./(nn1(ok).*nn2(ok)); cosab(~ok)=NaN;
nr(nn1>0)=nn2(nn1>0)./nn1(nn1>0); nr(nn1<=0)=NaN;
f3=figure('Color','w','Position',[100 100 1200 760]);
subplot(2,2,1); plotMap(Vert,Face,cosab,[-90 0]); cb=colorbar; ylabel(cb,'cos sim'); caxis([min(cosab) 1]); title('filter similarity (lat L)');
subplot(2,2,2); plotMap(Vert,Face,cosab,[90 0]); cb=colorbar; ylabel(cb,'cos sim'); caxis([min(cosab) 1]); title('filter similarity (lat R)');
ma=max(abs(log2(nr(isfinite(nr))))); cm=divmap(256);
subplot(2,2,3); plotMap(Vert,Face,log2(nr),[-90 0]); colormap(gca,cm); cb=colorbar; ylabel(cb,'log_2 norm ratio'); caxis([-ma ma]); title('gain redistribution (lat L)');
subplot(2,2,4); plotMap(Vert,Face,log2(nr),[90 0]); colormap(gca,cm); cb=colorbar; ylabel(cb,'log_2 norm ratio'); caxis([-ma ma]); title('gain redistribution (lat R)');
sgtitle('Per-vertex spatial-filter agreement: no-depth MNE vs Eigen-MNE');
saveFig(f3,fullfile(outDir,'kernel_agreement_maps'));

% ================= FIG 4: resolution (diagonal full + PSF subsample) =================
diag_nd=sum(K_nd.*G.',2); diag_eig=sum(K_eig.*G.',2);   % [nV x 1] cheap, no full R
nSub=min(3000,nV); rng_idx=randperm(nV,nSub);
Rnd_s=K_nd(rng_idx,:)*G; Reg_s=K_eig(rng_idx,:)*G;       % [nSub x nV]
psf_nd=zeros(nSub,1); psf_eig=zeros(nSub,1);
for k=1:nSub
    a=abs(Rnd_s(k,:)); psf_nd(k)=sum(a>0.5*abs(Rnd_s(k,rng_idx(k))));
    bb=abs(Reg_s(k,:)); psf_eig(k)=sum(bb>0.5*abs(Reg_s(k,rng_idx(k))));
end
f4=figure('Color','w','Position',[100 100 1100 420]);
subplot(1,2,1); scatter(diag_nd,diag_eig,4,[.5 .5 .5],'filled','MarkerFaceAlpha',0.3); hold on;
ld=[min([diag_nd;diag_eig]) max([diag_nd;diag_eig])]; plot(ld,ld,'k--'); grid on;
xlabel('MNE noDepth R_{ii}'); ylabel('Eigen R_{ii}'); cc2=corrcoef(diag_nd,diag_eig); title(sprintf('resolution diagonal r=%.3f',cc2(1,2)));
subplot(1,2,2); scatter(psf_nd,psf_eig,5,[.5 .5 .5],'filled','MarkerFaceAlpha',0.3); hold on;
lp=[0 max([psf_nd;psf_eig])]; plot(lp,lp,'k--'); grid on; axis equal tight;
xlabel('MNE noDepth PSF (#verts)'); ylabel('Eigen PSF (#verts)'); title(sprintf('PSF width (median %.0f vs %.0f)',median(psf_nd),median(psf_eig)));
sgtitle(sprintf('Resolution comparison (%d-vertex subsample)',nSub));
saveFig(f4,fullfile(outDir,'kernel_resolution.png'));

% ----- summary file -----
fid=fopen(fullfile(outDir,'kernel_comparison_summary.txt'),'w');
fprintf(fid,'Kernel comparison (TutorialAuditory, ico5)  sources=%d channels=%d\n',nV,nCh);
fprintf(fid,'SVD eff.rank99: MNEdepth=%d  MNEnoDepth=%d  Eigen=%d\n', r99(sv_std), r99(sv_nd), r99(sv_eig));
fprintf(fid,'no-depth vs Eigen: weight r=%.4f\n', cc(1,2));
fprintf(fid,'cos sim (noDepth vs Eigen): mean=%.4f median=%.4f min=%.4f  (cos<0.9: %.1f%%)\n', ...
    mean(cosab,'omitnan'), median(cosab,'omitnan'), min(cosab), 100*sum(cosab<0.9)/nV);
fprintf(fid,'norm ratio (Eigen/noDepth): median=%.4f\n', median(nr,'omitnan'));
fprintf(fid,'resolution diagonal r=%.4f ; PSF median noDepth=%.0f Eigen=%.0f\n', cc2(1,2), median(psf_nd), median(psf_eig));
fprintf(fid,'row-norm CV: MNEdepth=%.3f MNEnoDepth=%.3f\n', std(sqrt(sum(K_std.^2,2)))/mean(sqrt(sum(K_std.^2,2))), std(nn1)/mean(nn1));
fclose(fid);

% ----- assertions -----
assert(abs(cc(1,2))<=1 && cc(1,2)>0, 'no-depth vs Eigen weight correlation must be positive');
assert(median(cosab,'omitnan')>0.5, 'median filter cosine similarity unexpectedly low');
assert(r99(sv_eig)>0 && r99(sv_std)>0, 'SVD ranks must be positive');
assert(exist(fullfile(outDir,'kernel_svd_spectrum.png'),'file')==2,'svd figure not saved');
fprintf('SUMMARY: r(noDepth,Eigen)=%.3f  cosMed=%.3f  rank99 dep/nodep/eig=%d/%d/%d\n', ...
    cc(1,2), median(cosab,'omitnan'), r99(sv_std), r99(sv_nd), r99(sv_eig));
disp('ALL TESTS PASSED');
end

% ===== helpers (toolbox-free) =====
function v=getfield_d(s,f,d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end

function plotMap(V,F,data,va)
patch('Faces',F,'Vertices',V,'FaceVertexCData',data,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
material dull; camlight headlight; axis equal off; view(va); colormap(gca,parula(256));
end

function cmap=divmap(n)
half=floor(n/2); blue=[0.20 0.40 0.80]; red=[0.80 0.20 0.20]; white=[1 1 1];
c1=[linspace(blue(1),white(1),half)' linspace(blue(2),white(2),half)' linspace(blue(3),white(3),half)'];
c2=[linspace(white(1),red(1),n-half)' linspace(white(2),red(2),n-half)' linspace(white(3),red(3),n-half)'];
cmap=[c1;c2];
end

function saveFig(h,p)
if ~endsWith(p,'.png'); p=[p '.png']; end
print(h,p,'-dpng','-r150'); close(h);
end
