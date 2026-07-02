function tests = test_dirac_connectome_axes
tests = functiontests(localfunctions);
end
function test_axes_and_rowmap(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));   % controller launches a Dirac-dSPM session
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac-Connectome');
    axs = bst_eigen('Axes', struct('SurfaceFile',ax.SurfaceFile,'Variant','LB-Connectome','nModes',60,'TimeWindow',[0 .04],'SampleRate',100));
    nV = numel(axs.GlobalVertices{1});  K = size(axs.Phi{1},2);
    verifyEqual(t, size(ax.Phi{1}), [4*nV 3*K]);
    verifyEqual(t, numel(ax.Lambda{1}), 3*K);
    [C,kind] = bst_eigenfilter('Fiber', ax);
    verifyEqual(t, C, 4);  verifyEqual(t, kind, 'quaternion');
    % RowMap maps a 3-vector source into the quaternion imag slots for this variant
    F = zeros(3*nV,1);  [srcRows, dstRows, nrows, msg] = bst_eigenfilter('RowMap', F, ax, 1);
    verifyEmpty(t, msg);  verifyEqual(t, nrows, 4*nV);
end
