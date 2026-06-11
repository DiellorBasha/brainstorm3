function test_source_vector_glyphs()
% TEST_SOURCE_VECTOR_GLYPHS  Headless regression for the source-vector quiver helpers.
% Requires Brainstorm running (figure_3d on path). Tests pure subfunctions via dispatch.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % ===== ComputeSourceVectorGlyphs =====
    P   = [0 0 0; 1 0 0; 2 0 0];           % anchors
    Nrm = [0 0 1; 0 0 1; 0 0 1];           % +z normals
    V3  = [3 0 0;  0 0 0;  0 4 0];         % vtx1 +x(mag3), vtx2 zero, vtx3 +y(mag4)
    G = figure_3d('ComputeSourceVectorGlyphs', P, Nrm, V3, [], 0.5, 0.1);
    [nPass,nFail] = chk('glyph offset lifts base along +z', all(abs(G.Z - 0.1) < 1e-12), nPass,nFail);
    [nPass,nFail] = chk('glyph unit-normalized +x * scale', abs(G.U(1)-0.5)<1e-12 && abs(G.V(1))<1e-12 && abs(G.W(1))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('glyph zero vector -> zero arrow',   abs(G.U(2))<1e-12 && abs(G.V(2))<1e-12 && abs(G.W(2))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('glyph unit-normalized +y * scale',  abs(G.V(3)-0.5)<1e-12 && abs(G.U(3))<1e-12 && abs(G.W(3))<1e-12, nPass,nFail);
    G2 = figure_3d('ComputeSourceVectorGlyphs', P, Nrm, V3, [1 3], 1, 0);
    [nPass,nFail] = chk('glyph idx subselect count==2', numel(G2.X)==2, nPass,nFail);
    % Degenerate zero-length normal -> no lift (base stays at the anchor)
    Gz = figure_3d('ComputeSourceVectorGlyphs', [5 6 7], [0 0 0], [1 0 0], [], 1, 0.1);
    [nPass,nFail] = chk('glyph zero normal -> no offset', abs(Gz.X-5)<1e-12 && abs(Gz.Y-6)<1e-12 && abs(Gz.Z-7)<1e-12, nPass,nFail);

    fprintf('\n==== test_source_vector_glyphs: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_source_vector_glyphs: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
