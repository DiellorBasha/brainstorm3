function R = proto_connectome_spread(C, Sc, csvFile, accumCol, seedLabels)
% PROTO_CONNECTOME_SPREAD: test whether a PET tracer's longitudinal accumulation follows the
% structural connectome (template tractography). Two per-subject regressions on clean UCBERKELEY ADNI
% ROI SUVR (Desikan; the anatomical mismatch cancels in the longitudinal difference, so no per-subject
% anatomy is needed), each controlling for the region's own baseline:
%   (a) NEIGHBOUR : Delta ~ own_baseline + (Crw*baseline)         -- connectome-neighbour spread
%   (b) SEED      : Delta ~ own_baseline + connection-to-seed     -- spread from a Braak-origin seed
% Subjects are split at the median baseline accumCol into EARLY vs ADVANCED (staged spread).
%
% USAGE: R = proto_connectome_spread(C, Sc, csvFile, 'META_TEMPORAL_SUVR', {'entorhinal L','entorhinal R'})
%   C,Sc from proto_combined_laplacian (region x region connectome + Desikan scouts, same order).
%
% FINDINGS (28Jun2026 UCBERKELEY, 68-region Desikan connectome from a 105k-streamline ICBM152 template):
%   POOLED neighbour beta ~ 0 for both tracers (null) -- BUT staging/seed-focus reveal real spread:
%   (a) TAU neighbour beta: early -0.027, ADVANCED +0.030 (p=0.026) -> connectome spread emerges in the
%       advanced/spreading phase, diluted to null when pooled with early (anchoring-at-origin) subjects.
%   (b) SEED (entorhinal): TAU +0.022(early)/+0.050(adv) vs AMYLOID -0.134 -> tau accumulates in
%       entorhinal-CONNECTED regions (Braak circuit); amyloid ANTI-correlates (loads the neocortical DMN
%       far from the seed). tau-vs-amyloid gap ~0.18 >> SE(~0.01): robust, Braak-predicted direction.
%   Effects are weak in absolute terms (template connectome, Desikan, 6mm) but the tau/amyloid contrast
%   is strong. Answer: tau spreads via fibres from the entorhinal origin; amyloid does NOT (different
%   spatial driver). See [[face-method-rigor]]-style caveats: subject DWI + finer parcellation would sharpen it.
%
% Author: Diellor Basha, 2026 (prototype; template connectome)
    if nargin<5||isempty(seedLabels), seedLabels={'entorhinal L','entorhinal R'}; end
    T=readtable(csvFile,'VariableNamingRule','preserve'); vn=T.Properties.VariableNames;
    nReg=numel(Sc); Crw=C./max(sum(C,2),eps);
    col=strings(nReg,1);                                   % Desikan scout -> UCBERKELEY column
    for r=1:nReg, p=strsplit(Sc(r).Label); h="RH"; if strcmp(p{end},'L'), h="LH"; end
        col(r)=sprintf('CTX_%s_%s_SUVR',h,upper(p{1})); end
    have=ismember(col,vn); idx=find(have); cc=col(idx);
    Cs=Crw(idx,idx); Cs=Cs./max(sum(Cs,2),eps);
    seedReg=find(ismember({Sc.Label},seedLabels)); conn2seed=sum(C(idx,seedReg),2); seedIn=ismember(idx,seedReg);
    zc=@(x)(x-mean(x))/max(std(x),eps);
    subs=unique(string(T.PTID)); rec=[];   % [baseline-accum, neighbour-beta, seed-beta]
    for s=1:numel(subs)
        Ts=sortrows(T(strcmp(string(T.PTID),subs(s)),:),'SCANDATE');
        if height(Ts)<2 || (Ts.(accumCol)(end)-Ts.(accumCol)(1))<0.05, continue; end
        base=Ts{1,cellstr(cc)}'; last=Ts{end,cellstr(cc)}'; ok=base>0&last>0&isfinite(base)&isfinite(last);
        if sum(ok)<20, continue; end
        d=last-base;
        bN=[zc(base(ok)) zc(Cs(ok,:)*base)]\zc(d(ok));                  % neighbour
        ok2=ok&~seedIn; bS=[zc(base(ok2)) zc(conn2seed(ok2))]\zc(d(ok2)); % seed
        rec=[rec; Ts.(accumCol)(1) bN(2) bS(2)]; %#ok<AGROW>
    end
    thr=median(rec(:,1)); early=rec(:,1)<=thr; adv=rec(:,1)>thr;
    R=struct('nSubj',size(rec,1),'baseline',rec(:,1),'betaNeighbour',rec(:,2),'betaSeed',rec(:,3), ...
             'thr',thr,'idxEarly',early,'idxAdv',adv);
    fprintf('%s (%d accumulating subjects):\n', accumCol, R.nSubj);
    fprintf('  neighbour beta: early %+.3f | advanced %+.3f\n', mean(rec(early,2)), mean(rec(adv,2)));
    fprintf('  seed beta     : early %+.3f | advanced %+.3f  (seed=%s)\n', mean(rec(early,3)), mean(rec(adv,3)), strjoin(seedLabels,'+'));
end
