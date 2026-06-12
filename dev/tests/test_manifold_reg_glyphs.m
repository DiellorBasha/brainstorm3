function test_manifold_reg_glyphs()
% TEST_MANIFOLD_REG_GLYPHS  Headless regression for view_manifold_registration's
% pure RegSingGlyphs (sphere-pole lollipop geometry).
% Requires Brainstorm on path so view_manifold_registration('RegSingGlyphs', ...) dispatches.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % Two offset unit spheres. L centered at [0 1 0], R at [0 -1 0].
    % il = [1 2 3 4] (L), ir = [5 6 7 8] (R). Each: center + 4 unit directions.
    sphV = [ 0 1  1;  1 1  0; -1 1  0;  0 1 -1; ...    % L (center [0 1 0])
             0 -1 1;  1 -1 0; -1 -1 0;  0 -1 -1];      % R (center [0 -1 0])
    il = [1;2;3;4]; ir = [5;6;7;8];
    poleIdx = [1; 8];   % L-top [0 1 1], R-bottom [0 -1 -1]

    [Base, Tip] = view_manifold_registration('RegSingGlyphs', sphV, poleIdx, ir, il);
    [nPass,nFail] = chk('Base = sphV(poleIdx)', isequal(Base, sphV(poleIdx,:)), nPass,nFail);
    [nPass,nFail] = chk('count matches poleIdx', size(Tip,1)==2, nPass,nFail);
    % outward from each hemisphere center: L center [0 1 0], R center [0 -1 0]
    radL = sphV(1,:) - [0 1 0];   % [0 0 1]
    radR = sphV(8,:) - [0 -1 0];  % [0 0 -1]
    [nPass,nFail] = chk('L tip outward', dot(Tip(1,:)-Base(1,:), radL) > 0, nPass,nFail);
    [nPass,nFail] = chk('R tip outward', dot(Tip(2,:)-Base(2,:), radR) > 0, nPass,nFail);
    % the two pins point in opposite global z
    [nPass,nFail] = chk('L/R pins opposite z', sign(Tip(1,3)-Base(1,3)) == -sign(Tip(2,3)-Base(2,3)), nPass,nFail);
    % lift magnitude = 0.15 * radius (=1 here)
    [nPass,nFail] = chk('lift = 0.15*radius', abs(norm(Tip(1,:)-Base(1,:)) - 0.15) < 1e-12, nPass,nFail);

    fprintf('\n==== test_manifold_reg_glyphs: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_manifold_reg_glyphs: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
