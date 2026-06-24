function test_geodesic_move()
% TEST_GEODESIC_MOVE: the Area tool is gone from Scouts; the dynamics tool drives the scroll.
%
% USAGE:  test_geodesic_move   % Brainstorm running
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % T1: panel_scout no longer exposes the Area tool verbs (removed), Line tool stays
    haveArea = ~isempty(which('panel_scout'));
    areaGone = true;
    try, panel_scout('IsAreaToolActive');  areaGone = false; catch, areaGone = true; end %#ok<CTCH>
    lineStays = true;
    try, panel_scout('IsGeodesicToolActive'); catch, lineStays = false; end %#ok<CTCH>
    ok1 = haveArea && areaGone && lineStays;
    fprintf('T1 scout surgery: areaVerbGone=%d lineToolStays=%d => %s\n', areaGone, lineStays, PF{ok1+1});
    pass = pass && ok1;

    % T2: figure_3d source no longer references the removed scout Area verbs, and references the dynamics tool
    src = fileread(which('figure_3d'));
    noScoutArea = isempty(strfind(src, 'AreaToolScroll')) && isempty(strfind(src, 'IsAreaToolActive')); %#ok<STREMP>
    hasDynTool  = ~isempty(strfind(src, 'bst_geodesic_tool')); %#ok<STREMP>
    ok2 = noScoutArea && hasDynTool;
    fprintf('T2 figure_3d rewire: noScoutAreaRefs=%d hasDynTool=%d => %s\n', noScoutArea, hasDynTool, PF{ok2+1});
    pass = pass && ok2;

    % T3: with the tool OFF, OnScroll passes through (returns false) so the wheel still zooms
    bst_geodesic_tool('Toggle', 0);                 % ensure inactive
    passthrough = bst_geodesic_tool('OnScroll', -1);
    ok3 = islogical(passthrough) && (passthrough == false);
    fprintf('T3 scroll passthrough (tool off): consumed=%d => %s\n', passthrough, PF{ok3+1});
    pass = pass && ok3;

    % T4: scout pick-mode clears the dynamics pick flag (bidirectional mutual exclusion)
    src2 = fileread(which('panel_scout'));
    ok4 = ~isempty(strfind(src2, 'isDynamicsGeodesicPick')); %#ok<STREMP>
    fprintf('T4 mutual-exclusion: scout clears dyn flag=%d => %s\n', ok4, PF{ok4+1});
    pass = pass && ok4;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
