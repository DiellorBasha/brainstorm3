function tests = test_scalogram_modes
tests = functiontests(localfunctions);
end
function test_scalogram_from_mode_c(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig,'DynamicsOverlay');
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac');
    nV = 0; for h=1:2, nV=max(nV,max(ax.GlobalVertices{h}(:))); end
    [C, gvAll] = panel_bst_dynamics('i_apply_projection', st, ax, D, 1:8, nV);
    verifyEqual(t, numel(C), 2);
    verifyEqual(t, size(C{1},1), size(ax.Lambda{1},1));   % coefficients = inverse mode count, per hemi
    verifyEqual(t, size(C{1},2), 8);
end
