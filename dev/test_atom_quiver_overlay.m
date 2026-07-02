function tests = test_atom_quiver_overlay
tests = functiontests(localfunctions);
end

% NOTE: SetAtomField runs its real scalar-paint body (bst_colormaps + i_dynamics_overlay) before the
% quiver block. On a bare test figure we (1) seed the 'Colormap' appdata AddColormapToFigure needs, and
% (2) pass W=[] so the atom paint (i_dynamics_overlay) early-returns on empty AtomField -- isolating the
% V3/QuiverVectorOverride contract this test targets without a full rendered surface.
function i_scaffold(hFig)
    setappdata(hFig, 'Colormap', struct('AllTypes',{{}}, 'Type','', 'DisplayUnits',[]));
end

function test_v3_sets_override(tc)
    hFig = figure('Visible','off');  c = onCleanup(@() close(hFig)); %#ok<NASGU>
    setappdata(hFig, 'DynamicsOverlay', struct('Op','none','iTess',1,'nV',4,'srcDS',[],'srcResult',[], ...
        'AtomField',[],'AtomGV',[],'AtomSigned',false,'AtomWin',[]));
    setappdata(hFig, 'Surface', struct('hPatch',[],'DataSource',struct('Type','Source','FileName','x'), ...
        'ShowSourceVectors',false));
    i_scaffold(hFig);
    V3 = [0 0 1; 0 0 1; 0 0 0; 0 0 0];
    view_dynamics('SetAtomField', hFig, [], (1:4)', false, V3);
    verifyEqual(tc, getappdata(hFig,'QuiverVectorOverride'), V3);
end

function test_empty_v3_clears(tc)
    hFig = figure('Visible','off');  c = onCleanup(@() close(hFig)); %#ok<NASGU>
    setappdata(hFig, 'DynamicsOverlay', struct('Op','none','iTess',1,'nV',4,'srcDS',[],'srcResult',[], ...
        'AtomField',[],'AtomGV',[],'AtomSigned',false,'AtomWin',[]));
    setappdata(hFig, 'Surface', struct('hPatch',[],'DataSource',struct('Type','Source','FileName','x'), ...
        'ShowSourceVectors',true));
    i_scaffold(hFig);
    setappdata(hFig, 'QuiverVectorOverride', [0 0 1]);
    view_dynamics('SetAtomField', hFig, [], (1:4)', false, []);
    verifyEmpty(tc, getappdata(hFig,'QuiverVectorOverride'));
end
