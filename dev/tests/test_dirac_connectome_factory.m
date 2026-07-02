function test_dirac_connectome_factory()
% Factory Dirac-Connectome eigenbasis == the reference lift of the LB-Connectome eigenbasis.
% Requires a live Brainstorm session with the 'preventad' protocol loaded (run via MATLAB-MCP
% or `matlab -batch`); the controller runs this, it is not executed as part of static checks.
if ~brainstorm('status'), brainstorm nogui; end
iP = bst_get('Protocol','preventad'); gui_brainstorm('SetCurrentProtocol', iP);
surf = 'sub-MTL0005/tess_cortex_pial_low.mat';

% Reference: build the LB-Connectome basis directly, then lift it by hand (the pre-refactor
% panel recipe).
axc = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','LB-Connectome','nModes',60));
[Pq, Lq, Mq] = bst_lift_connectome_dirac(axc.Phi{1}, axc.Lambda{1}(:), axc.Mass{1});

% Factory: tess_operators/tess_eigen build+lift 'Dirac-Connectome' internally; the panel now
% only calls bst_eigen('Axes','Dirac-Connectome') (no inline bst_lift_connectome_dirac call).
axd = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac-Connectome','nModes',60));

assert(isequal(size(axd.Phi{1}), size(Pq)) && max(abs(axd.Phi{1}(:)-Pq(:)))<1e-10, 'Phi mismatch');
assert(isequal(size(axd.Lambda{1}), size(Lq)) && max(abs(axd.Lambda{1}(:)-Lq(:)))<1e-10, 'Lambda mismatch');
assert(isequal(size(axd.Mass{1}), size(Mq)) && max(abs(axd.Mass{1}(:)-Mq(:)))<1e-10, 'Mass mismatch');

fs = bst_eigen('FieldSpec', axd);
assert(strcmp(fs.field_type,'quaternion'), 'field_type');

fprintf('PASS\n');
end
