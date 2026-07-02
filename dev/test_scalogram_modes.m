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
    % Scalogram must actually RUN on the quaternion Dirac coefficients (regression: it assumed a scalar
    % basis and errored on W(gv,:,m)=Wm because manifold_ift returns the [4nV] quaternion field).
    lmax = max([ax.Lambda{1}(:); ax.Lambda{2}(:)]);  N = 4;
    gCell = cell(1,N); for m=1:N, gCell{m} = bst_eigfilter_design_itersine(struct('member',m,'Nf',N,'lmax',lmax)); end
    scal = bst_eigenwavelet('Scalogram', ax, gCell, C);
    verifyEqual(t, size(scal.W,3), N);                    % per-band reconstructed fields
    verifyEqual(t, size(scal.W,1), nV);                   % PER-VERTEX magnitude (not 4*nV quaternion rows)
    verifyEqual(t, size(scal.energy,3), N);
    verifyGreaterThan(t, max(scal.centers), 0.5*sqrt(lmax)); % spans the full inverse mode range, not a 60-mode cap
end
