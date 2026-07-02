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
    % amplitude recon == dSPM in DIRECTION *per vertex* (they differ only by the per-vertex SIR scale, so a
    % global flattened cosine would be dragged below 1 by that reweighting -- must check per-vertex).
    Jgt = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0)); % [3nV x nWin]
    V1 = V3(:,:,1);                                % [nV x 3] amplitude recon, frame 1
    G1 = reshape(Jgt(:,1),3,[])';                  % [nV x 3] dSPM, frame 1
    cosv = zeros(nV,1);
    for v=1:nV, a=V1(v,:); b=G1(v,:); na=norm(a); nb=norm(b); if na>0&&nb>0, cosv(v)=dot(a,b)/(na*nb); end, end
    verifyGreaterThan(t, median(cosv(cosv~=0)), 0.999);   % per-vertex direction identical to dSPM
end
