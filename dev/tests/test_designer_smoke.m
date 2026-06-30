% test_designer_smoke - the designer's realiser path end-to-end on a real cortex (no GUI)
SURF = 'sub-MTL0002/tess_cortex_pial_low.mat';
ax = bst_eigen('Axes', struct('SurfaceFile',SURF, 'Variant','Laplace-Beltrami', 'nModes',50, ...
                              'TimeWindow',[0 0.99], 'SampleRate',100));
seed = ax.GlobalVertices{1}(1);
% diffusion (ts) atom
[Wd,gv] = bst_eigenfilter('Atom', ax, 'diffusion', struct('lmax',max(ax.Lambda{1}),'tau',0.3), seed);
assert(isequal(size(Wd),[numel(gv) ax.nT]) && all(isfinite(Wd(:))), 'ts (diffusion) atom end-to-end');
% heat (static) atom
[Ws,~] = bst_eigenfilter('Atom', ax, 'heat', struct('lmax',max(ax.Lambda{1}),'t',0.2), seed);
assert(max(abs(Ws(:,1)-Ws(:,end))) < 1e-9, 'static (heat) atom constant in time');
disp('OK');
