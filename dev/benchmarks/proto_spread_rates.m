function R = proto_spread_rates(csvFile, accumCol, hotspotLabels, Sc, col, SurfaceFile)
% PROTO_SPREAD_RATES: physical-unit rates of tracer accumulation/spread from longitudinal UCBERKELEY
% ADNI ROI SUVR + scan dates. Three quantities, each with real units:
%   accumulation rate (SUVR/yr) = mean(last-base)/dt           -- how fast intensity grows
%   growth rate r (1/yr)        = mean(log(last/base))/dt       -- exponential rate -> doubling log(2)/r
%   spatial spread velocity (mm/yr) = d/dt of the pathology-weighted spatial spread-radius
%                                     radius = sqrt(sum w|pos-ctr|^2/sum w), w=max(SUVR-1,0)
% The first two are AMPLITUDE; the third is EUCLIDEAN spatial expansion. Network (connectome) spread
% is NOT captured here (fiber-connected regions are Euclidean-distant) - that is a redistribution rate
% (connectome-predicted SUVR/yr, see proto_connectome_spread), not a mm/yr front velocity.
%
% USAGE: R = proto_spread_rates(csv,'META_TEMPORAL_SUVR',{'entorhinal L','entorhinal R'},Sc,col,wf)
%   Sc,col = Desikan scouts + UCBERKELEY-column mapping (from proto_connectome_spread setup); wf = a
%   manifold cortex for region centroids (e.g. sub-MTL0002 white).
%
% FINDINGS (28Jun2026 UCBERKELEY, n=423 amyloid / 182 tau accumulating):
%   AMYLOID: 0.026 SUVR/yr cortical-mean (0.034 precuneus), r=0.020/yr (doubling 34yr), +0.16 mm/yr
%   TAU:     0.013 SUVR/yr cortical-mean (0.030 entorhinal=focal), r=0.010/yr (doubling 67yr), +0.05 mm/yr
%   => spread is AMPLITUDE-dominated (SUVR/yr), not a fast Euclidean front (mm/yr tiny for both); amyloid
%   broad+slowly-expanding, tau focal+seed-anchored (its spread is network, not Euclidean).
%
% Author: Diellor Basha, 2026 (prototype)
    T=readtable(csvFile,'VariableNamingRule','preserve'); vn=T.Properties.VariableNames;
    have=ismember(col,vn); idx=find(have); cc=col(idx);
    keyCol=col(ismember({Sc.Label},hotspotLabels)); keyCol=keyCol(ismember(keyCol,vn));
    sW=in_tess_bst(SurfaceFile); ScM=sW.Atlas(find(strcmp({sW.Atlas.Name},'Desikan-Killiany'),1)).Scouts;
    cent=zeros(numel(idx),3);
    for k=1:numel(idx), j=find(strcmp({ScM.Label},Sc(idx(k)).Label),1); if ~isempty(j), cent(k,:)=mean(sW.Vertices(ScM(j).Vertices,:),1)*1000; end; end
    subs=unique(string(T.PTID)); rate=[]; rg=[]; krate=[]; vel=[];
    for s=1:numel(subs)
        Ts=sortrows(T(strcmp(string(T.PTID),subs(s)),:),'SCANDATE');
        if height(Ts)<2 || (Ts.(accumCol)(end)-Ts.(accumCol)(1))<0.05, continue; end
        dt=years(Ts.SCANDATE(end)-Ts.SCANDATE(1)); if dt<0.5, continue; end
        base=Ts{1,cellstr(cc)}'; last=Ts{end,cellstr(cc)}'; ok=base>0&last>0&isfinite(base)&isfinite(last);
        rate(end+1)=mean(last(ok)-base(ok))/dt; rg(end+1)=mean(log(last(ok)./base(ok)))/dt; %#ok<AGROW>
        krate(end+1)=(mean(Ts{end,cellstr(keyCol)})-mean(Ts{1,cellstr(keyCol)}))/dt; %#ok<AGROW>
        rr=nan(1,2); XX={base,last};
        for q=1:2, w=max(XX{q}-1,0); w(~isfinite(w))=0; if sum(w)>eps, ctr=(w'*cent)/sum(w); rr(q)=sqrt(sum(w.*sum((cent-ctr).^2,2))/sum(w)); end; end
        vel(end+1)=(rr(2)-rr(1))/dt; %#ok<AGROW>
    end
    md=@(x)median(x,'omitnan');
    R=struct('accumRate',rate,'growthRate',rg,'hotspotRate',krate,'spreadVel',vel, ...
             'accumRate_med',md(rate),'doublingYr',log(2)/md(rg),'hotspotRate_med',md(krate),'spreadVel_med',md(vel));
    fprintf('%s (n=%d): accum %.4f SUVR/yr | r %.4f/yr (doubling %.0fyr) | hotspot %.4f SUVR/yr | spread %+.2f mm/yr\n', ...
        accumCol, numel(rate), R.accumRate_med, md(rg), R.doublingYr, R.hotspotRate_med, R.spreadVel_med);
end
