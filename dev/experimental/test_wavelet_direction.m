function test_wavelet_direction()
% Pure test of the local-frame direction embedding d = cos(t)(cos(p)U+sin(p)V)+sin(t)N.
% Authors: Diellor Basha, 2026
    nFail = 0;
    U = [1 0 0]; V = [0 1 0]; N = [0 0 1];     % orthonormal frame
    % theta=+90 -> exactly the (outward) normal
    d = panel_wavelet_designer('EmbedDirection', 0, 90, U, V, N);
    nFail = nFail + chk('tilt +90 = normal', max(abs(d - N)) < 1e-12);
    % theta=0, phi=0 -> U
    d = panel_wavelet_designer('EmbedDirection', 0, 0, U, V, N);
    nFail = nFail + chk('tilt 0, phi 0 = U', max(abs(d - U)) < 1e-12);
    % theta=0, phi=90 -> V
    d = panel_wavelet_designer('EmbedDirection', 90, 0, U, V, N);
    nFail = nFail + chk('tilt 0, phi 90 = V', max(abs(d - V)) < 1e-12);
    % general: unit length, correct formula
    p = 37; t = 20;
    d = panel_wavelet_designer('EmbedDirection', p, t, U, V, N);
    expect = cosd(t)*(cosd(p)*U + sind(p)*V) + sind(t)*N; expect = expect/norm(expect);
    nFail = nFail + chk('general formula', max(abs(d - expect)) < 1e-12);
    nFail = nFail + chk('unit length', abs(norm(d)-1) < 1e-12);
    fprintf('\n==== test_wavelet_direction: %d failed ====\n', nFail);
    if nFail > 0, error('test_wavelet_direction FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
