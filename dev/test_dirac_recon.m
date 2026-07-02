function tests = test_dirac_recon
tests = functiontests(localfunctions);
end
function test_recon_matches_amplitude_direction(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig,'DynamicsOverlay');
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac');
    iWin = 1:4;
    [cCell,~] = panel_bst_dynamics('i_mode_coeffs', st, D, iWin);
    [V3, mag] = panel_bst_dynamics('i_dirac_recon', ax, cCell);
    nV = size(V3,1);  verifyEqual(t, size(V3,3), 4);  verifyEqual(t, size(mag), [nV 4]);
    % direction equals GetResultsValues (dSPM) per vertex at frame 1
    Jgt = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0)); % [3nV x nWin]
    v1 = reshape(V3(:,:,1)',[],1);                 % [3nV x 1] recon frame 1
    g1 = Jgt(:,1);
    cang = dot(v1,g1)/(norm(v1)*norm(g1)+eps);
    verifyGreaterThan(t, cang, 0.999);             % same overall direction (amplitude vs dSPM = per-vtx scale)
end
