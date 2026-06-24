function test_bst_geodesic_tool()
% TEST_BST_GEODESIC_TOOL: the dynamics heat-disk tool (seed/grow/state + draw/clear).
%
% USAGE:  test_bst_geodesic_tool   % Brainstorm running, TutorialAuditory loaded
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % cortex surface for the geometry
    sSubject = bst_get('Subject', 'Subject01');
    if isempty(sSubject) || isempty(sSubject.iCortex)
        fprintf('SKIPPED (no Subject01 cortex)\n');  fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    Surf = in_tess_bst(SurfaceFile, 0);
    vi = round(size(Surf.Vertices,1)/3);

    % ---------- T1: Seed -> GetState ----------
    bst_geodesic_tool('Seed', SurfaceFile, vi);
    st = bst_geodesic_tool('GetState');
    n0 = 0;  if ~isempty(st), n0 = numel(st.vertices); end
    ok1 = ~isempty(st) && (st.seed==vi) && (n0>0) && isequal(size(st.pos),[1 3]) ...
       && strcmp(st.SurfaceFile, SurfaceFile) && isequal(st.pos, Surf.Vertices(vi,:));
    fprintf('T1 seed: verts=%d seed=%d pos=%d => %s\n', n0, (st.seed==vi), isequal(st.pos,Surf.Vertices(vi,:)), PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2: Grow (up = -1 grows; down = +N shrinks; radius floors) ----------
    bst_geodesic_tool('Grow', -1);  s2 = bst_geodesic_tool('GetState');  nUp = numel(s2.vertices);
    bst_geodesic_tool('Grow', +10); s3 = bst_geodesic_tool('GetState');  nDn = numel(s3.vertices);
    ok2 = (nUp >= n0) && (nDn <= nUp) && (s3.radius >= 0.003 - 1e-12) && (s2.radius > s3.radius);
    fprintf('T2 grow: n0=%d up=%d down=%d radiusFloor=%d => %s\n', n0, nUp, nDn, (s3.radius>=0.003-1e-12), PF{ok2+1});
    pass = pass && ok2;

    % ---------- T3: Draw/Clear overlay on a figure ----------
    hFig = view_surface(SurfaceFile);  drawnow;
    bst_geodesic_tool('Seed', SurfaceFile, vi);
    bst_geodesic_tool('Draw', hFig);   drawnow;
    nPatch = numel(findobj(hFig, 'Tag', 'GeodesicToolDisk'));
    bst_geodesic_tool('Clear', hFig);  drawnow;
    nAfter = numel(findobj(hFig, 'Tag', 'GeodesicToolDisk'));
    ok3 = (nPatch==1) && (nAfter==0);
    fprintf('T3 draw/clear: patch=%d cleared=%d => %s\n', nPatch, (nAfter==0), PF{ok3+1});
    pass = pass && ok3;
    if ishandle(hFig), close(hFig); end

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
