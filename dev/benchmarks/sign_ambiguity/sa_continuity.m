function C = sa_continuity(F, SulciMap, ResultsFile, opts)
% SA_CONTINUITY: Component 2 of the sign-ambiguity continuity test. Applies an
% unconstrained MEG kernel, picks the peak sample, and compares phase
% continuity across facing sulcal-wall pairs under three frames.
%
% USAGE:  C = sa_continuity(F, SulciMap, ResultsFile, opts)
%   F           : struct from sa_frames.
%   SulciMap    : [nV x 1] binary.
%   ResultsFile : link file ('link|kernel|data') or full results file with nComponents=3.
%                 Bare kernel files (no paired data) are not supported (LoadFull=1 required).
%   opts        : MaxDist, NormalDot, Nring (pair detection); Window [t0 t1] s
%                 (peak search); MagRatio (continuity conditioning).
%
% OUTPUT C: J (peak current [nV x 3]), tPeak, pairs, conditioned mask, per-frame
% phase discontinuities (circular distance in [0,pi]) and the constrained
% sign-flip rate.

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'MaxDist'),  opts.MaxDist  = 0.003; end
if ~isfield(opts,'NormalDot'),opts.NormalDot= -0.7;  end
if ~isfield(opts,'Nring'),    opts.Nring    = 3;     end
if ~isfield(opts,'Window'),   opts.Window   = [0.05 0.15]; end
if ~isfield(opts,'MagRatio'), opts.MagRatio = 0.5;   end

% ===== Apply kernel and pick the peak sample =====
Results = in_bst_results(ResultsFile, 1);   % LoadFull=1 -> Kernel*Data
assert(isequal(Results.nComponents, 3), 'sa_continuity requires nComponents=3.');
IGA  = Results.ImageGridAmp;                % [3nV x nTime]
Time = Results.Time(:)';
nV   = size(F.Nv, 1);
assert(size(IGA,1) == 3*nV, ...
    'sa_continuity:sizeMismatch', ...
    'IGA rows (%d) do not match 3*nV (%d); surface/results mismatch.', size(IGA,1), 3*nV);
inWin = Time >= opts.Window(1) & Time <= opts.Window(2);
if ~any(inWin), inWin = true(size(Time)); end
gfp = sqrt(sum(IGA.^2, 1));                 % global field power over sources
gfp(~inWin) = -inf;
[~, ti] = max(gfp);
J = reshape(IGA(:, ti), 3, []).';           % [nV x 3] component-interleaved
C.tPeak = Time(ti);

% ===== Sulcal-wall pairs =====
pairs = sa_sulcal_walls(F.Vtx, F.Nv, SulciMap, F.VertConn, opts);
i = pairs(:,1); j = pairs(:,2);

% ===== Conditioning on physical continuity =====
Ji = J(i,:); Jj = J(j,:);
mi = sqrt(sum(Ji.^2,2)); mj = sqrt(sum(Jj.^2,2));
cosang = sum(Ji.*Jj,2) ./ max(mi.*mj, eps);
magOk  = min(mi,mj) ./ max(max(mi,mj), eps) > opts.MagRatio;
keep   = (cosang > 0) & magOk;
i = i(keep); j = j(keep);
C.pairs = [i j];
C.nPairs = numel(i);
if C.nPairs == 0
    C.signFlipRate = NaN;
    C.phaseDisc = struct('fiedler_median',NaN,'tang_median',NaN,'xyz_median',NaN);
    C.phaseDiscRate = struct('fiedler',NaN,'tang',NaN,'xyz',NaN);
    C.df = []; C.dt = []; C.dx = [];
    C.J = J; return;
end

% ===== Constrained sign-flip rate =====
s = sum(J .* F.Nv, 2);                       % constrained scalar J . n
C.signFlipRate = mean(sign(s(i)) ~= sign(s(j)));

% ===== Per-frame phase discontinuity (circular distance) =====
phaseOf = @(fr) atan2(sum(J .* fr.e2, 2), sum(J .* fr.e1, 2));
circDist = @(a,b) abs(atan2(sin(a-b), cos(a-b)));
pf = phaseOf(F.fiedler); pt = phaseOf(F.tang); px = phaseOf(F.xyz);
fOk = sqrt(sum(F.fiedler.e1(i,:).^2,2)) > 1e-9 & ...
      sqrt(sum(F.fiedler.e1(j,:).^2,2)) > 1e-9;
df = circDist(pf(i), pf(j));
df(~fOk) = NaN;   % exclude singularity vertices from Fiedler stats
dt = circDist(pt(i), pt(j));
dx = circDist(px(i), px(j));
C.phaseDisc = struct('fiedler_median', median(df(fOk)), ...
                     'tang_median',    median(dt), ...
                     'xyz_median',     median(dx));
% A pair is "phase-discontinuous" if the circular distance exceeds pi/2.
C.phaseDiscRate = struct('fiedler', mean(df(fOk) > pi/2), ...
                         'tang',    mean(dt > pi/2), ...
                         'xyz',     mean(dx > pi/2));
% Carry per-pair vectors for figures.
C.df = df; C.dt = dt; C.dx = dx; C.J = J;
end
