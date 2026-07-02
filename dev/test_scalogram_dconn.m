function tests = test_scalogram_dconn
% Live test (needs a booted Dirac-dSPM Dynamics session): the Dirac-Connectome projection returns
% vector (3K) coefficients and bst_eigenwavelet('Scalogram') reduces them to a per-vertex magnitude.
tests = functiontests(localfunctions);
end
function test_scalogram(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig, 'DynamicsOverlay');
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac-Connectome');
    nV = numel(ax.GlobalVertices{1});
    [C, ~] = panel_bst_dynamics('i_apply_projection', st, ax, D, 1:8, nV);
    verifyEqual(t, size(C{1},1), size(ax.Lambda{1},1));   % 3K vector coeffs
    lmax = max(ax.Lambda{1}(:)); N = 4;
    gC = cell(1,N); for m=1:N, gC{m}=bst_eigfilter_design_itersine(struct('member',m,'Nf',N,'lmax',lmax)); end
    scal = bst_eigenwavelet('Scalogram', ax, gC, C);
    verifyEqual(t, size(scal.W,1), nV);                  % per-vertex magnitude (not 4nV)
    verifyEqual(t, size(scal.energy,3), N);
end
