function test_eigfilter_pure
% Verify the bst_eigfilter library: factories return handles, vector scale -> bank,
% 'meta' specs, registry list/info/dispatch/error, evaluate (handle + bank), compose.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

lam = [0; 1; 4; 9; 16; 25];

% ---- heat factory: handle, analytic value, t->0 ~ 1 ----
g = bst_eigfilter_design_heat(struct('t', 0.1));
assert(isa(g, 'function_handle'), 'factory must return a handle.');
assert(max(abs(g(lam) - exp(-0.1*lam))) < 1e-12, 'heat value wrong.');
g0 = bst_eigfilter_design_heat(struct('t', 1e-12));
assert(max(abs(g0(lam) - 1)) < 1e-6, 'heat t->0 must be ~1.');

% ---- vector scale -> cell bank ----
gb = bst_eigfilter_design_heat(struct('t', [0.1 0.5 1.0]));
assert(iscell(gb) && numel(gb) == 3, 'vector scale must return a 3-cell bank.');
assert(max(abs(gb{2}(lam) - exp(-0.5*lam))) < 1e-12, 'bank element wrong.');

% ---- meta ----
m = bst_eigfilter_design_heat('meta');
assert(strcmp(m.name,'heat') && isfield(m,'params') && isfield(m,'priorAdmissible'), 'meta malformed.');

% ---- evaluate: handle -> [K x 1], bank -> [K x Nf] ----
h1 = bst_eigfilter_evaluate(g, lam);
assert(isequal(size(h1), [numel(lam) 1]), 'evaluate handle shape.');
hB = bst_eigfilter_evaluate(gb, lam);
assert(isequal(size(hB), [numel(lam) 3]), 'evaluate bank shape.');

% ---- compose: product of two handles ----
g2 = bst_eigfilter_design_heat(struct('t', 0.2));
gc = bst_eigfilter_compose(g, g2);
assert(max(abs(gc(lam) - g(lam).*g2(lam))) < 1e-12, 'compose must be pointwise product.');

% ---- registry: dispatch, list, info, unknown error ----
gr = bst_eigfilter_kernel('heat', struct('t', 0.1));
assert(max(abs(gr(lam) - exp(-0.1*lam))) < 1e-12, 'registry dispatch wrong.');
names = bst_eigfilter_kernel('list');
assert(any(strcmp(names, 'heat')), 'list must include heat.');
mi = bst_eigfilter_kernel('info', 'heat');
assert(strcmp(mi.name, 'heat'), 'info must return heat meta.');
threw = false;
try, bst_eigfilter_kernel('no_such_kernel'); catch, threw = true; end
assert(threw, 'unknown kernel must error.');

lp = [1; 2; 4; 8; 16];   % strictly positive for power/log
% flat
gf = bst_eigfilter_design_flat();
assert(all(bst_eigfilter_evaluate(gf, lp) == 1), 'flat must be ones.');
% power: lambda^-alpha, decreasing
gp = bst_eigfilter_design_power(struct('alpha',1));
assert(max(abs(gp(lp) - lp.^(-1))) < 1e-12, 'power value wrong.');
% log: -log(lambda), needs lambda in (0,1) for positivity (prior scales it)
gl = bst_eigfilter_design_log();
assert(max(abs(gl([0.1;0.5]) - (-log([0.1;0.5])))) < 1e-12, 'log value wrong.');
% inverse_heat: clamped
gih = bst_eigfilter_design_inverse_heat(struct('t',0.1,'maxgain',5));
assert(max(gih(lp)) <= 5 + 1e-12, 'inverse_heat must clamp at maxgain.');
% tikhonov
gt = bst_eigfilter_design_tikhonov(struct('beta',2));
assert(max(abs(gt(lp) - 1./(1+2*lp))) < 1e-12, 'tikhonov value wrong.');
% ideal band [2 8] inclusive
gi = bst_eigfilter_design_ideal(struct('band',[2 8]));
assert(isequal(gi(lp), double(lp>=2 & lp<=8)), 'ideal mask wrong.');
% registry sees all of them
nm = bst_eigfilter_kernel('list');
for k = {'flat','power','log','inverse_heat','tikhonov','ideal'}
    assert(any(strcmp(nm, k{1})), sprintf('list missing %s.', k{1}));
end

lpp = [1; 2; 4; 8; 16];
% matern: (kappa^2 + l)^-nu, decreasing, positive
gm = bst_eigfilter_design_matern(struct('kappa',1,'nu',1.5));
assert(max(abs(gm(lpp) - (1 + lpp).^(-1.5))) < 1e-12, 'matern value wrong.');
assert(all(diff(gm(lpp)) < 0), 'matern must be decreasing.');
% mexhat band-pass: zero at 0, peak interior, decay; vector t -> bank
gh = bst_eigfilter_design_mexhat(struct('t',0.1));
assert(abs(gh(0)) < 1e-12, 'mexhat must be 0 at lambda=0.');
ll = (0:0.01:50)'; v = gh(ll);
assert(v(1) < max(v) && v(end) < max(v), 'mexhat must peak in the interior.');
ghb = bst_eigfilter_design_mexhat(struct('t',[0.05 0.1 0.2]));
assert(iscell(ghb) && numel(ghb)==3, 'mexhat vector t must return a bank.');
mh = bst_eigfilter_design_mexhat('meta');
assert(mh.priorAdmissible == false, 'mexhat must be flagged not prior-admissible.');
% diffgauss band-pass: non-negative for t1<t2, zero at 0
gd = bst_eigfilter_design_diffgauss(struct('t1',0.1,'t2',0.4));
assert(abs(gd(0)) < 1e-12, 'diffgauss must be 0 at lambda=0.');
assert(all(gd(lpp) >= -1e-12), 'diffgauss must be non-negative for t1<t2.');
md = bst_eigfilter_design_diffgauss('meta');
assert(md.priorAdmissible == false, 'diffgauss must be flagged not prior-admissible.');

% ---- meta contract: every registered kernel exposes the full field set ----
% (mechanically enforces the contract the future filter-design GUI relies on)
reqFields = {'name','display','params','bandpass','priorAdmissible'};
allNames  = bst_eigfilter_kernel('list');
assert(numel(allNames) == 10, 'expected 10 registered kernels.');
for k = allNames
    mt = bst_eigfilter_kernel('info', k{1});
    for f = reqFields
        assert(isfield(mt, f{1}), sprintf('meta for %s missing field %s.', k{1}, f{1}));
    end
    assert(strcmp(mt.name, k{1}), sprintf('meta.name mismatch for %s.', k{1}));
    assert(islogical(mt.bandpass) || isnumeric(mt.bandpass), 'meta.bandpass must be logical/numeric.');
    assert(islogical(mt.priorAdmissible) || isnumeric(mt.priorAdmissible), 'meta.priorAdmissible must be logical/numeric.');
end
% only mexhat and diffgauss are analysis-only (priorAdmissible=false)
inadmissible = {};
for k = allNames
    mt = bst_eigfilter_kernel('info', k{1});
    if ~mt.priorAdmissible; inadmissible{end+1} = k{1}; end %#ok<AGROW>
end
assert(isequal(sort(inadmissible), {'diffgauss','mexhat'}), 'only diffgauss and mexhat must be prior-inadmissible.');

disp('ALL TESTS PASSED');
end
