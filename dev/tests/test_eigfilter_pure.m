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

disp('ALL TESTS PASSED');
end
