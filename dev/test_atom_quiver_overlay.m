function tests = test_atom_quiver_overlay
tests = functiontests(localfunctions);
end

function test_v3_sets_override(tc)
    hFig = figure('Visible','off');  c = onCleanup(@() close(hFig));
    setappdata(hFig, 'DynamicsOverlay', struct('Op','none','iTess',1,'nV',4,'srcDS',[],'srcResult',[], ...
        'AtomField',[],'AtomGV',[],'AtomSigned',false,'AtomWin',[]));
    setappdata(hFig, 'Surface', struct('hPatch',[],'DataSource',struct('Type','Source','FileName','x'), ...
        'ShowSourceVectors',false));
    V3 = [0 0 1; 0 0 1; 0 0 0; 0 0 0];
    view_dynamics('SetAtomField', hFig, [1;1;0;0], (1:4)', false, V3);
    verifyEqual(tc, getappdata(hFig,'QuiverVectorOverride'), V3);
end

function test_empty_v3_clears(tc)
    hFig = figure('Visible','off');  c = onCleanup(@() close(hFig));
    setappdata(hFig, 'DynamicsOverlay', struct('Op','none','iTess',1,'nV',4,'srcDS',[],'srcResult',[], ...
        'AtomField',[],'AtomGV',[],'AtomSigned',false,'AtomWin',[]));
    setappdata(hFig, 'Surface', struct('hPatch',[],'DataSource',struct('Type','Source','FileName','x'), ...
        'ShowSourceVectors',true));
    setappdata(hFig, 'QuiverVectorOverride', [0 0 1]);
    view_dynamics('SetAtomField', hFig, [1;1;0;0], (1:4)', false, []);
    verifyEmpty(tc, getappdata(hFig,'QuiverVectorOverride'));
end
