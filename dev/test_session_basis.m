function tests = test_session_basis
tests = functiontests(localfunctions);
end
function test_uses_inverse_basis(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac');
    D  = getappdata(st.hFig,'DynamicsOverlay');
    src = panel_bst_dynamics('i_src_resultfile', D);
    R = in_bst_results(src,0,'Eigenvalues','ModeHemisphere','DiracEigenFile');
    verifyEqual(t, size(ax.Lambda{1},1), sum(R.ModeHemisphere==1));   % inverse's mode count, not 60
    verifyGreaterThan(t, size(ax.Lambda{1},1), 60);
    verifyEqual(t, size(ax.Phi{1},1), 4*numel(ax.GlobalVertices{1})); % quaternion layout
    verifyTrue(t, issorted(ax.Lambda{1}));
end
