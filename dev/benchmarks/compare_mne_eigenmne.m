function out = compare_mne_eigenmne(protocolName)
% COMPARE_MNE_EIGENMNE: Vertex-by-vertex comparison of standard MNE vs Eigen-MNE
% source maps on the folder-averaged deviant/standard evoked responses, plus a
% direct comparison of the two imaging kernels. Produces figures + a summary.
%
% USAGE:  out = compare_mne_eigenmne('TutorialAuditory')
if nargin<1 || isempty(protocolName); protocolName='TutorialAuditory'; end
iProt=bst_get('Protocol',protocolName); gui_brainstorm('SetCurrentProtocol',iProt);
outDir=fullfile(fileparts(mfilename('fullpath')),'mne_vs_eigenmne'); if ~exist(outDir,'dir'); mkdir(outDir); end

% ---- locate a study holding BOTH kernels + BOTH folder averages ----
sStudies=bst_get('ProtocolStudies'); T=[];
for iS=1:numel(sStudies.Study)
  s=sStudies.Study(iS); if isempty(s.Result)||isempty(s.Data); continue; end
  rc=string({s.Result.Comment}); dc=string({s.Data.Comment});
  hasK = any(rc=="MN: MEG(Constr) 2018") && any(rc=="Eigen-MNE: MEG 2018");
  hasA = any(startsWith(dc,"Avg: deviant")) && any(startsWith(dc,"Avg: standard"));
  if hasK && hasA; T=iS; break; end
end
if isempty(T); error('No study with both kernels and both folder averages. Run the folder average first.'); end
s=sStudies.Study(T); rc=string({s.Result.Comment}); dc=string({s.Data.Comment});
iMN=find(rc=="MN: MEG(Constr) 2018"); iMN=iMN(end);
iEG=find(rc=="Eigen-MNE: MEG 2018"); iEG=iEG(end);
R1=in_bst_results(s.Result(iMN).FileName,0); R2=in_bst_results(s.Result(iEG).FileName,0);
K1=R1.ImagingKernel; K2=R2.ImagingKernel;          % [nV x nCh]
V=in_tess_bst(R1.SurfaceFile);

conds={'deviant','standard'};
out=struct(); out.protocol=protocolName;
rT=struct();  % per-time correlation per condition
for ic=1:2
  iD=find(startsWith(dc,"Avg: "+conds{ic})); iD=iD(1);
  D=in_bst_data(s.Data(iD).FileName); Time=D.Time;
  J1=K1*D.F(R1.GoodChannel,:); J2=K2*D.F(R2.GoodChannel,:);   % [nV x nT]
  % per-time signed Pearson correlation across vertices
  A=J1-mean(J1,1); B=J2-mean(J2,1);
  r = sum(A.*B,1) ./ sqrt(sum(A.^2,1).*sum(B.^2,1));
  gfp1=sqrt(mean(J1.^2,1)); gfp2=sqrt(mean(J2.^2,1));
  rT.(conds{ic})=struct('Time',Time,'r',r,'gfp1',gfp1,'gfp2',gfp2,'J1',J1,'J2',J2);
end

% ===== F1: per-time correlation + GFP, both conditions =====
f1=figure('Color','w','Position',[100 100 1100 480]);
for ic=1:2
  C=rT.(conds{ic});
  subplot(2,2,ic); plot(C.Time, C.r,'LineWidth',1.2); ylim([-0.2 1]); grid on;
  ylabel('vertex corr r(t)'); title(sprintf('%s : MNE vs Eigen-MNE', conds{ic})); yline(0,'k:');
  subplot(2,2,ic+2); plot(C.Time, C.gfp1,'b', C.Time, C.gfp2,'r','LineWidth',1); grid on;
  xlabel('Time (s)'); ylabel('GFP'); legend('MNE','Eigen-MNE','Location','best');
end
sgtitle('Vertex-wise agreement over time (signed source values)');
print(f1,fullfile(outDir,'f1_corr_timecourse.png'),'-dpng','-r150'); close(f1);

% ===== F2: vertex scatter at the auditory evoked peak (standard, cleaner) =====
condDisp='standard'; C=rT.(condDisp);
win = C.Time>=0.05 & C.Time<=0.25; gfpw=C.gfp1; gfpw(~win)=-inf; [~,tpk]=max(gfpw);
x=C.J1(:,tpk); y=C.J2(:,tpk);
f2=figure('Color','w','Position',[100 100 560 520]);
lim=max(abs([x;y])); nb=120; edges=linspace(-lim,lim,nb);
H=histcounts2(x,y,edges,edges); imagesc(edges,edges,log10(H'+1)); axis xy equal tight; hold on;
plot([-lim lim],[-lim lim],'w--'); colormap(hot); colorbar;
xlabel('MNE source value'); ylabel('Eigen-MNE source value');
title(sprintf('%s @ t=%.3fs  (r=%.3f)', condDisp, C.Time(tpk), C.r(tpk)));
print(f2,fullfile(outDir,'f2_vertex_scatter_peak.png'),'-dpng','-r150'); close(f2);

% ===== F3: cortex maps at deviant peak (MNE | Eigen-MNE | difference) =====
m1=C.J1(:,tpk); m2=C.J2(:,tpk); md=m1-m2;
f3=figure('Color','w','Position',[100 100 1320 420]); maps={m1,m2,md}; tit={'MNE','Eigen-MNE','MNE - Eigen-MNE'};
for k=1:3
  subplot(1,3,k); cl=max(abs(maps{k}));
  patch('Vertices',V.Vertices,'Faces',V.Faces,'FaceVertexCData',maps{k},'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
  axis equal off; view(-90,90); camlight headlight; caxis([-cl cl]); colormap(jet); title(tit{k});
end
sgtitle(sprintf('%s evoked source map @ t=%.3fs', condDisp, C.Time(tpk)));
print(f3,fullfile(outDir,'f3_cortex_maps_peak.png'),'-dpng','-r150'); close(f3);

% ===== F4: direct ImagingKernel comparison =====
% per-vertex (row) correlation of the 272-channel mapping
A=K1-mean(K1,2); B=K2-mean(K2,2);
rowr = sum(A.*B,2) ./ sqrt(sum(A.^2,2).*sum(B.^2,2));
n1=sqrt(sum(K1.^2,2)); n2=sqrt(sum(K2.^2,2));        % per-vertex row norm
ccAll=corrcoef(K1(:),K2(:)); kAll=ccAll(1,2);
relDiff=norm(K1-K2,'fro')/norm(K1,'fro');
f4=figure('Color','w','Position',[100 100 1100 440]);
subplot(1,2,1); histogram(rowr,50); grid on; xlabel('per-vertex kernel-row correlation'); ylabel('# vertices');
title(sprintf('row corr: median=%.3f', median(rowr)));
subplot(1,2,2); loglog(n1,n2,'.','MarkerSize',3); hold on; ll=[min([n1;n2]) max([n1;n2])]; loglog(ll,ll,'k--'); axis equal tight; grid on;
xlabel('||MNE row||'); ylabel('||Eigen-MNE row||'); title('per-vertex kernel gain');
sgtitle(sprintf('ImagingKernel comparison  (overall corr=%.3f, rel.Fro.diff=%.3f)', kAll, relDiff));
print(f4,fullfile(outDir,'f4_kernel_comparison.png'),'-dpng','-r150'); close(f4);

% ---- summary ----
out.peakTime=C.Time(tpk); out.peakCorr=C.r(tpk);
out.medianRowCorr=median(rowr); out.kernelOverallCorr=kAll; out.kernelRelFroDiff=relDiff;
out.devCorrRange=[min(rT.deviant.r) max(rT.deviant.r)];
out.figDir=outDir;
fprintf('SUMMARY:\n');
fprintf('  evoked (standard) peak t=%.3fs  vertex r=%.3f\n', out.peakTime, out.peakCorr);
fprintf('  kernel overall corr=%.3f  rel Fro diff=%.3f  median row corr=%.3f\n', kAll, relDiff, median(rowr));
fprintf('  figures -> %s\n', outDir);
end
