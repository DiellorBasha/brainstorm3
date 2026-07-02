function tests = test_mode_coeffs
tests = functiontests(localfunctions);
end
function tc = i_launch(~)
    % helper: launch a Dirac-dSPM Dynamics session, return st + D
end
function test_shape_and_split(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st), 'launch a Dirac dSPM Dynamics session first');
    D  = getappdata(st.hFig,'DynamicsOverlay');
    assertTrue(t, panel_bst_dynamics('i_is_dirac_dspm', D));
    iWin = 1:5;
    [cCell, meta] = panel_bst_dynamics('i_mode_coeffs', st, D, iWin);
    verifyEqual(t, numel(cCell), 2);
    verifyEqual(t, size(cCell{1},2), 5);
    % per-hemi mode counts match ModeHemisphere
    src = panel_bst_dynamics('i_src_resultfile', D);
    R = in_bst_results(src,0,'ModeHemisphere');
    verifyEqual(t, size(cCell{1},1), sum(R.ModeHemisphere==1));
    verifyEqual(t, size(cCell{2},1), sum(R.ModeHemisphere==2));
    % each hemi ascending lambda
    verifyTrue(t, issorted(meta.Eigenvalues{1}));
end
