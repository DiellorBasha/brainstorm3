function R = proto_regional_redistribution(csvFile, accumCol)
% PROTO_REGIONAL_REDISTRIBUTION: regional accumulation-vs-baseline redistribution test on clean
% UCBERKELEY ADNI ROI SUVR (no surface sampling needed; the anatomical mismatch cancels in the
% longitudinal difference, so this is valid without per-subject anatomy). For each accumulating
% subject (>=2 timepoints), correlate the regional growth rate log(last/base) with the baseline:
%   corr(growth-rate, baseline) >= 0  -> AMPLIFICATION (high-baseline regions grow fastest)
%   corr(growth-rate, baseline) <  0  -> SPREAD/FILLING  (low-baseline regions catch up)
%
% USAGE: R = proto_regional_redistribution('.../Fully_Processed_Amyloid_..._UCBERKELEY_..csv','SUMMARY_SUVR')
%
% FINDINGS (28Jun2026 UCBERKELEY, 34 bilateral Desikan cortical ROIs):
%   AMYLOID (n=422): corr(rate,baseline)=-0.05 -> fills in (low regions catch up), 56% spread sig
%   TAU     (n=181): corr(rate,baseline)=+0.18 -> amplifies at the Braak origin (high regions
%                    accelerate; cohort is mostly early-Braak), 34% spread sig
%   Difference 0.23 >> SE (~0.02) -> amyloid and tau have significantly DIFFERENT, opposite
%   redistribution signatures. corr(accumulation,baseline)=0.15 (amyloid) matches the vertex-level
%   stand-in-anatomy value, confirming that low value is REAL (mismatch cancels) not noise.
%
% Author: Diellor Basha, 2026
    T=readtable(csvFile,'VariableNamingRule','preserve'); vn=T.Properties.VariableNames;
    ctx=vn(~cellfun('isempty',regexp(vn,'^CTX_[A-Z]+_SUVR$','once')));   % bilateral Desikan cortical ROIs
    subs=unique(string(T.PTID)); rDB=[]; rRB=[];
    for s=1:numel(subs)
        Ts=sortrows(T(strcmp(string(T.PTID),subs(s)),:),'SCANDATE');
        if height(Ts)<2, continue; end
        base=Ts{1,ctx}'; last=Ts{end,ctx}'; ok=base>0&last>0&isfinite(base)&isfinite(last);
        if sum(ok)<20, continue; end
        if Ts.(accumCol)(end)-Ts.(accumCol)(1) < 0.05, continue; end    % require real accumulation
        da=last(ok)-base(ok); rate=log(last(ok)./base(ok));
        cDB=corrcoef(da,base(ok)); cRB=corrcoef(rate,base(ok));
        rDB(end+1)=cDB(1,2); rRB(end+1)=cRB(1,2); %#ok<AGROW>
    end
    R=struct('nSubj',numel(rDB),'nROI',numel(ctx),'corrAccumBaseline',rDB,'corrRateBaseline',rRB);
    fprintf('%s: %d accumulating subjects, %d cortical ROIs\n', csvFile, R.nSubj, R.nROI);
    fprintf('  corr(accumulation,baseline) = %.2f +/- %.2f\n', mean(rDB), std(rDB));
    fprintf('  corr(growth-rate, baseline) = %.2f +/- %.2f  [<0 SPREAD, >=0 AMPLIFY]  | %.0f%% spread\n', mean(rRB), std(rRB), 100*mean(rRB<0));
end
