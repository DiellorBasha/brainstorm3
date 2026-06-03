function St = bench_stats(csvPath, outDir)
% BENCH_STATS: Aggregate benchmark rows -> median/IQR per (method,regime,snr,K);
% paired Wilcoxon signed-rank (each eig method vs each standard method, matched on
% anatomy/regime/snr/replicate at the eig method's largest K). Writes stats.csv + stats.md.
if nargin < 2 || isempty(outDir); outDir = fileparts(csvPath); end
T = readtable(csvPath);
T.method = string(T.method); T.regime = string(T.regime); T.anatomy = string(T.anatomy);

eigMethods = ["eig_mne_log","eig_dspm_log"];
stdMethods = ["wmne","dspm","sloreta"];
regimes = unique(T.regime,'stable');
Kmax = max(T.K(~isnan(T.K)));

methods = [stdMethods, eigMethods];
S = table();
for ir = 1:numel(regimes)
    for im = 1:numel(methods)
        m = methods(im);
        if ismember(m, eigMethods)
            sel = T.method==m & T.regime==regimes(ir) & T.K==Kmax;
        else
            sel = T.method==m & T.regime==regimes(ir);
        end
        v = T.locerror_mm(sel);
        if isempty(v); continue; end
        row = table(string(regimes(ir)), m, median(v), iqr_local(v), mean(v), numel(v), ...
            'VariableNames', {'regime','method','median_locerror_mm','iqr_locerror_mm','mean_locerror_mm','n'});
        S = [S; row]; %#ok<AGROW>
    end
end

Cmp = table();
for ir = 1:numel(regimes)
    for ie = 1:numel(eigMethods)
        eTab = T(T.method==eigMethods(ie) & T.regime==regimes(ir) & T.K==Kmax, :);
        for is = 1:numel(stdMethods)
            sTab = T(T.method==stdMethods(is) & T.regime==regimes(ir), :);
            [d, n] = paired_diffs(eTab, sTab, 'locerror_mm');
            if n < 1; continue; end
            p = local_signrank(d);
            row = table(string(regimes(ir)), eigMethods(ie), stdMethods(is), median(d), n, p, ...
                'VariableNames', {'regime','eig_method','std_method','median_diff_mm','n_pairs','p_value'});
            Cmp = [Cmp; row]; %#ok<AGROW>
        end
    end
end

Ksweep = table();
Kvals = unique(T.K(~isnan(T.K)));
for ir = 1:numel(regimes)
    for ie = 1:numel(eigMethods)
        for ik = 1:numel(Kvals)
            v = T.locerror_mm(T.method==eigMethods(ie) & T.regime==regimes(ir) & T.K==Kvals(ik));
            if isempty(v); continue; end
            Ksweep = [Ksweep; table(string(regimes(ir)), eigMethods(ie), Kvals(ik), median(v), ...
                'VariableNames',{'regime','eig_method','K','median_locerror_mm'})]; %#ok<AGROW>
        end
    end
end

St = struct('summary',S,'compare',Cmp,'ksweep',Ksweep);
writetable(S,   fullfile(outDir,'stats.csv'));
writetable(Cmp, fullfile(outDir,'stats_compare.csv'));
writetable(Ksweep, fullfile(outDir,'stats_ksweep.csv'));

fid = fopen(fullfile(outDir,'stats.md'),'w');
fprintf(fid,'# Eigenmode accuracy benchmark - statistics\n\n## Median LocError (mm) by method x regime\n\n');
fprintf(fid,'| regime | method | median | IQR | n |\n|---|---|---|---|---|\n');
for i=1:height(S)
    fprintf(fid,'| %s | %s | %.2f | %.2f | %d |\n', S.regime(i), S.method(i), ...
        S.median_locerror_mm(i), S.iqr_locerror_mm(i), S.n(i));
end
fprintf(fid,'\n## Paired Wilcoxon (eig vs standard), LocError\n\n');
fprintf(fid,'| regime | eig | vs | median diff (mm) | n | p |\n|---|---|---|---|---|---|\n');
for i=1:height(Cmp)
    fprintf(fid,'| %s | %s | %s | %+.2f | %d | %.4g |\n', Cmp.regime(i), Cmp.eig_method(i), ...
        Cmp.std_method(i), Cmp.median_diff_mm(i), Cmp.n_pairs(i), Cmp.p_value(i));
end
fclose(fid);
end

function q = iqr_local(v); q = quantile_local(v,0.75) - quantile_local(v,0.25); end

function q = quantile_local(v, p)
v = sort(v(:)); n = numel(v);
if n==1; q=v; return; end
h = (n-1)*p + 1; lo = floor(h); hi = min(lo+1,n);
q = v(lo) + (h-lo)*(v(hi)-v(lo));
end

function [d, n] = paired_diffs(eTab, sTab, valVar)
ke = strcat(string(eTab.anatomy),'|',string(eTab.regime),'|',string(eTab.snr_db),'|',string(eTab.replicate));
ks = strcat(string(sTab.anatomy),'|',string(sTab.regime),'|',string(sTab.snr_db),'|',string(sTab.replicate));
[tf, loc] = ismember(ke, ks);
d = eTab.(valVar)(tf) - sTab.(valVar)(loc(tf));
n = numel(d);
end

function p = local_signrank(d)
% Two-sided Wilcoxon signed-rank p-value via normal approximation with tie correction.
d = d(d~=0); n = numel(d);
if n < 1; p = NaN; return; end
R = local_tiedrank(abs(d));
Wp = sum(R(d>0));
mu = n*(n+1)/4;
[~,~,g] = unique(abs(d));
tieTerm = 0;
for k = 1:max(g); tk = sum(g==k); tieTerm = tieTerm + (tk^3 - tk); end
sigma = sqrt(n*(n+1)*(2*n+1)/24 - tieTerm/48);
if sigma == 0; p = 1; return; end
z = (Wp - mu) / sigma;
p = erfc(abs(z)/sqrt(2));
end

function r = local_tiedrank(x)
[xs, ord] = sort(x(:)); r = zeros(numel(x),1); i = 1; n = numel(xs);
while i <= n
    j = i; while j < n && xs(j+1)==xs(i); j = j+1; end
    r(ord(i:j)) = (i+j)/2; i = j+1;
end
end
