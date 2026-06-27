function S = summarize_pvc_3way(ResFile)
% SUMMARIZE_PVC_3WAY: aggregate batch_pvc_3way results across subjects/tracers.
% Reports per-pair cortex r/CCC (mean +/- SD, split by tracer) and a POOLED analysis
% (all cortical ROIs from all subjects, each pre-normalized to its own cerebellum).
%
% USAGE:  S = summarize_pvc_3way()  % default results file
%
% Author: Diellor Basha, 2026

    here = bst_fileparts(mfilename('fullpath'));
    if (nargin < 1) || isempty(ResFile), ResFile = fullfile(here,'pvc_3way_results.mat'); end
    L = load(ResFile,'Res'); Res = L.Res;
    nE = numel(Res);
    pairNames = {Res(1).metrics.pair};           % 1=OURS-MGX, 2=MGX-GTM, 3=OURS-GTM
    shortp = {'OURS-MGX (same method)','MGX-GTM (PETSurfer)','OURS-GTM (cross)'};

    % per-entry cortex r/ccc
    rC = zeros(nE,3); cC = zeros(nE,3); trac = cell(nE,1);
    for i=1:nE
        rC(i,:) = [Res(i).metrics.r_ctx];
        cC(i,:) = [Res(i).metrics.ccc_ctx];
        trac{i} = Res(i).tracer;
    end
    isAmy = strcmp(trac,'18FNAV4694'); isTau = strcmp(trac,'18Fflortaucipir');

    fprintf('\n================= 3-WAY PVC SUMMARY (n=%d subj-tracer) =================\n', nE);
    fprintf('%-26s %16s %16s %16s\n','pair','amyloid r (m+/-sd)','tau r (m+/-sd)','ALL r (m+/-sd)');
    for k=1:3
        fprintf('%-26s   %5.2f +/- %4.2f     %5.2f +/- %4.2f     %5.2f +/- %4.2f\n', shortp{k}, ...
            mean(rC(isAmy,k)),std(rC(isAmy,k)), mean(rC(isTau,k)),std(rC(isTau,k)), mean(rC(:,k)),std(rC(:,k)));
    end
    fprintf('%-26s %16s %16s %16s\n','pair (CCC)','amyloid','tau','ALL');
    for k=1:3
        fprintf('%-26s   %5.2f +/- %4.2f     %5.2f +/- %4.2f     %5.2f +/- %4.2f\n', shortp{k}, ...
            mean(cC(isAmy,k)),std(cC(isAmy,k)), mean(cC(isTau,k)),std(cC(isTau,k)), mean(cC(:,k)),std(cC(:,k)));
    end

    % POOLED cortical ROIs across all subjects
    PO=[]; PM=[]; PT=[]; Ptr={};
    for i=1:nE
        rows = Res(i).rows; ic = strcmp({rows.class},'cortex');
        PO=[PO,[rows(ic).ours]]; PM=[PM,[rows(ic).mgx]]; PT=[PT,[rows(ic).gtm]]; %#ok<AGROW>
        Ptr=[Ptr, repmat(Res(i).tracer,1,nnz(ic))]; %#ok<AGROW>
    end
    PO=PO(:); PM=PM(:); PT=PT(:); pAmy=strcmp(Ptr(:),'18FNAV4694');
    cr=@(a,b) local_corr(a,b); cc=@(a,b) local_ccc(a,b);
    fprintf('\nPOOLED cortical ROIs (N=%d): \n', numel(PO));
    fprintf('  OURS-MGX  r=%.3f CCC=%.3f | MGX-GTM r=%.3f CCC=%.3f | OURS-GTM r=%.3f CCC=%.3f\n', ...
        cr(PO,PM),cc(PO,PM), cr(PM,PT),cc(PM,PT), cr(PO,PT),cc(PO,PT));

    % outliers (cortex r<0.3 any pair)
    bad = find(any(rC<0.3,2));
    if ~isempty(bad)
        fprintf('\n[outliers cortex r<0.3]:\n');
        for j=bad(:)', fprintf('   %s/%s : OURS-MGX=%.2f MGX-GTM=%.2f OURS-GTM=%.2f\n', Res(j).subject, Res(j).tracer, rC(j,1),rC(j,2),rC(j,3)); end
    end

    S.rC=rC; S.cC=cC; S.trac=trac; S.pooled=struct('O',PO,'M',PM,'T',PT,'isAmy',pAmy);

    % ===== figure =====
    f=figure('Visible','off','Position',[60 60 1500 460]);
    % (1) per-subject r by pair x tracer
    subplot(1,3,1); hold on;
    x=1:3; w=0.18;
    mA=mean(rC(isAmy,:)); sA=std(rC(isAmy,:)); mT=mean(rC(isTau,:)); sT=std(rC(isTau,:));
    bar(x-w, mA, 0.32,'FaceColor',[.2 .4 .8]); bar(x+w, mT, 0.32,'FaceColor',[.85 .4 .2]);
    errorbar(x-w,mA,sA,'k.','LineWidth',1); errorbar(x+w,mT,sT,'k.','LineWidth',1);
    for k=1:3
        scatter((x(k)-w)*ones(nnz(isAmy),1)+0.03*randn(nnz(isAmy),1), rC(isAmy,k),18,[.1 .2 .5],'filled','MarkerFaceAlpha',.6);
        scatter((x(k)+w)*ones(nnz(isTau),1)+0.03*randn(nnz(isTau),1), rC(isTau,k),18,[.5 .2 .1],'filled','MarkerFaceAlpha',.6);
    end
    set(gca,'XTick',1:3,'XTickLabel',{'OURS-MGX','MGX-GTM','OURS-GTM'},'XTickLabelRotation',15);
    ylabel('cortex r'); ylim([-0.6 1]); yline(0,'k:'); grid on; legend({'amyloid','tau'},'Location','southwest');
    title('per-subject cortex r by pair');
    % (2) pooled OURS vs MGX
    subplot(1,3,2); local_scatter(PM,PO,pAmy,'MGX (PETSurfer MG)','OURS (PETPVE12 MG)','POOLED OURS vs MGX (same method)');
    % (3) pooled OURS vs GTM
    subplot(1,3,3); local_scatter(PT,PO,pAmy,'GTM (PETSurfer)','OURS (PETPVE12 MG)','POOLED OURS vs GTM');
    sgtitle(sprintf('3-way PVC across %d subjects x 2 tracers', nE/2),'Interpreter','none');
    print(f, fullfile(here,'summarize_pvc_3way.png'),'-dpng','-r110'); close(f);
    fprintf('\nFigure: %s\n', fullfile(here,'summarize_pvc_3way.png'));
end

function local_scatter(x,y,isAmy,xl,yl,ti)
    hold on;
    scatter(x(isAmy),y(isAmy),14,[.2 .4 .8],'filled','MarkerFaceAlpha',.5);
    scatter(x(~isAmy),y(~isAmy),14,[.85 .4 .2],'filled','MarkerFaceAlpha',.5);
    lim=[0 max([x;y])*1.05]; plot(lim,lim,'k--'); xlim(lim); ylim(lim); axis square; grid on;
    xlabel(xl); ylabel(yl);
    title(sprintf('%s\nr=%.2f CCC=%.2f', ti, local_corr(x,y), local_ccc(x,y)));
    legend({'amyloid','tau','identity'},'Location','northwest');
end

function c = local_corr(x,y), cc=corrcoef(x,y); c=cc(1,2); end
function c = local_ccc(x,y)
    mx=mean(x); my=mean(y); vx=var(x,1); vy=var(y,1); sxy=mean((x-mx).*(y-my));
    c = 2*sxy/(vx+vy+(mx-my)^2);
end
