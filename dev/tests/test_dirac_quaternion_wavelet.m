function test_dirac_quaternion_wavelet()
% TEST_DIRAC_QUATERNION_WAVELET: full-quaternion end-to-end for the Dirac wavelet.
%
% Proves the Dirac wavelet carries the FULL quaternion (real w retained) through
% Analysis -> Synthesis, and that atoms steer EXACTLY by right-quaternion multiplication:
%   - Analysis output is the full-quaternion field [4*nV x nT x M] (not the 3-vector);
%   - the real (w) component is genuinely nonzero (would be lost by the old 3-slot map);
%   - tight-frame (itersine) Synthesis(Analysis(.)) == the band-limited quaternion
%     projection, to machine precision;
%   - ToVec/ToQuat are consistent and ToQuat is pure (w=0);
%   - Steer reproduces RE-SEEDING: Steer(atom(q), q0) == atom(q (x) q0)  [the steering law];
%   - the scalar Laplace-Beltrami path is unchanged (real, native layout, exact RT).
%
% Auto-detects Dirac + Laplace-Beltrami eigen nodes under Subject01; SKIPs if absent.
%
% Authors: Diellor Basha, 2026

    nFail = 0;
    chk = @(name, cond) i_chk(name, cond);

    efD = i_find_variant('Dirac');
    efL = i_find_variant('Laplace-Beltrami');
    if isempty(efD)
        fprintf('SKIP test_dirac_quaternion_wavelet: no Dirac eigen node found.\n'); return;
    end
    fprintf('Dirac node: %s\n', efD);

    ED = in_bst_eigen(efD);  OD = in_bst_operator(ED.OperatorFile);
    nVg = max(cellfun(@(v) max([v(:);0]), ED.GlobalVertices));
    lams = ED.Lambda(~cellfun(@isempty, ED.Lambda));
    lrange = [min(cellfun(@min,lams)), max(cellfun(@max,lams))];
    frT = bst_eigenwavelet('Design', 'itersine', 6, lrange);    % tight frame
    M   = numel(frT.g);
    nT  = 3;

    % physical 3-vector input (w=0), random
    F3 = randn(3*nVg, nT);

    % ---- (1) Analysis output is the FULL quaternion field ----
    W = bst_eigenwavelet('Analysis', F3, ED, OD, frT);
    nFail = nFail + chk('Analysis output is [4*nV x nT x M]', isequal(size(W), [4*nVg, nT, M]));
    nFail = nFail + chk('W is real', isreal(W));

    % ---- (2) the real (w) component is genuinely carried ----
    Wq = reshape(W, 4, nVg, nT*M);
    wmax = max(abs(Wq(1,:,:)), [], 'all');  vmax = max(abs(Wq(2:4,:,:)), [], 'all');
    nFail = nFail + chk('w component nonzero (carried end-to-end)', wmax > 1e-6 * vmax);

    % ---- (3) tight-frame round-trip == band-limited quaternion projection ----
    Frec = bst_eigenwavelet('Synthesis', W, ED, OD, frT);
    nFail = nFail + chk('Synthesis output is [4*nV x nT]', isequal(size(Frec), [4*nVg, nT]));
    proj = i_bandlimit_quat(F3, ED, OD);                         % direct ΦΦ*B projection (w=0 input)
    rt = norm(Frec(:) - proj(:)) / max(norm(proj(:)), eps);
    nFail = nFail + chk('tight-frame RT == projection (machine exact)', rt < 1e-9);

    % ---- (4) ToVec / ToQuat ----
    V = bst_eigenwavelet('ToVec', W, ED);
    nFail = nFail + chk('ToVec shape [3*nV x nT x M]', isequal(size(V), [3*nVg, nT, M]));
    Q = bst_eigenwavelet('ToQuat', F3, ED);
    Qq = reshape(Q, 4, nVg, nT);
    nFail = nFail + chk('ToQuat is pure (w=0)', max(abs(Qq(1,:,:)),[],'all') == 0);
    nFail = nFail + chk('ToVec(ToQuat(V))==V', isequal(bst_eigenwavelet('ToVec', Q, ED), F3));

    % ---- (5) Steer applies right-quaternion multiplication per quaternion ----
    q0 = [cos(pi/5); 0.3; -0.4; 0.2];  q0 = q0/norm(q0);
    Ws = bst_eigenwavelet('Steer', W, q0, ED);
    Wsq = reshape(Ws, 4, nVg, nT*M);  ref = zeros(size(Wsq));
    for k=1:nVg, for t=1:nT*M, ref(:,k,t) = i_qmul(Wq(:,k,t), q0); end, end %#ok<ALIGN>
    nFail = nFail + chk('Steer == per-quaternion right-mult', norm(Wsq(:)-ref(:))/max(norm(Wq(:)),eps) < 1e-12);

    % ---- (6) STEERING LAW: Steer(atom(q),q0) == atom(q (x) q0) ----
    seed = ED.GlobalVertices{1}(round(numel(ED.GlobalVertices{1})/2));
    v = [1;0;0];  qs = [0; v];                                   % pure-imaginary seed
    A_pure = bst_eigenwavelet('Atom', ED, OD, frT, seed, v);     % API atom (3-vec seed)
    A_man  = i_manual_atom(ED, OD, frT, seed, qs);               % manual ref, same pure seed
    nFail = nFail + chk('Atom API == manual (pure seed)', norm(A_pure(:)-A_man(:))/max(norm(A_man(:)),eps) < 1e-10);
    A_steer = bst_eigenwavelet('Steer', A_pure, q0, ED);
    A_reseed = i_manual_atom(ED, OD, frT, seed, i_qmul(qs, q0));  % atom for the rotated seed
    nFail = nFail + chk('STEER law: Steer(atom(q),q0)==atom(q*q0)', ...
                        norm(A_steer(:)-A_reseed(:))/max(norm(A_reseed(:)),eps) < 1e-10);

    % ---- (7) scalar Laplace-Beltrami path unchanged ----
    if ~isempty(efL)
        EL = in_bst_eigen(efL);  OL = in_bst_operator(EL.OperatorFile);
        nVl = max(cellfun(@(v_) max([v_(:);0]), EL.GlobalVertices));
        lamsL = EL.Lambda(~cellfun(@isempty, EL.Lambda));
        frL = bst_eigenwavelet('Design', 'itersine', 6, [min(cellfun(@min,lamsL)) max(cellfun(@max,lamsL))]);
        FL = randn(nVl, 2);
        WL = bst_eigenwavelet('Analysis', FL, EL, OL, frL);
        nFail = nFail + chk('LBO Analysis native layout [nV x nT x M]', isequal(size(WL),[nVl,2,numel(frL.g)]) && isreal(WL));
        RL = bst_eigenwavelet('Synthesis', WL, EL, OL, frL);
        pL = i_bandlimit_scalar(FL, EL, OL);
        nFail = nFail + chk('LBO tight-frame RT == projection', norm(RL(:)-pL(:))/max(norm(pL(:)),eps) < 1e-9);
    end

    fprintf('\n==== test_dirac_quaternion_wavelet: %d failed ====\n', nFail);
    if nFail > 0; error('test_dirac_quaternion_wavelet FAILED'); end
    disp('ALL TESTS PASSED');
end

% ===== helpers =====
function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end

function r = i_qmul(a, b)
    r = [ a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4); ...
          a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3); ...
          a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2); ...
          a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1) ];
end

function P = i_bandlimit_quat(F3, EM, OM)
    % full-quaternion band-limited projection of the pure (w=0) embedding of F3
    nVg = max(cellfun(@(v) max([v(:);0]), EM.GlobalVertices));
    P = zeros(4*nVg, size(F3,2));
    for h = 1:numel(EM.Phi)
        Phi = EM.Phi{h}; if isempty(Phi); continue; end
        idx = EM.GlobalVertices{h}(:); n = numel(idx);
        gIn  = reshape([(idx-1)*3+1,(idx-1)*3+2,(idx-1)*3+3].',[],1);
        lIn  = reshape([(0:n-1)*4+2;(0:n-1)*4+3;(0:n-1)*4+4],[],1);
        gOut = reshape([(idx-1)*4+1,(idx-1)*4+2,(idx-1)*4+3,(idx-1)*4+4].',[],1);
        U = zeros(4*n, size(F3,2)); U(lIn,:) = F3(gIn,:);
        Up = Phi*(Phi'*(OM.Mass{h}*U));
        P(gOut,:) = Up;
    end
end

function P = i_bandlimit_scalar(F, EM, OM)
    P = zeros(size(F));
    for h = 1:numel(EM.Phi)
        Phi = EM.Phi{h}; if isempty(Phi); continue; end
        gv = EM.GlobalVertices{h}(:);
        P(gv,:) = Phi*(Phi'*(OM.Mass{h}*F(gv,:)));
    end
end

function A = i_manual_atom(EM, OM, frame, seedGlobal, q4)
    % full-quaternion atom for a quaternion seed q4 at global vertex seedGlobal
    nVg = max(cellfun(@(v) max([v(:);0]), EM.GlobalVertices));
    M = numel(frame.g);
    A = zeros(4*nVg, 1, M);
    for h = 1:numel(EM.Phi)
        Phi = EM.Phi{h}; if isempty(Phi); continue; end
        idx = EM.GlobalVertices{h}(:); n = numel(idx);
        loc = find(idx == seedGlobal, 1);
        gOut = reshape([(idx-1)*4+1,(idx-1)*4+2,(idx-1)*4+3,(idx-1)*4+4].',[],1);
        U = zeros(4*n, 1);
        if ~isempty(loc); U((loc-1)*4 + (1:4)) = q4; end
        H = bst_eigenwavelet('Evaluate', frame, EM.Lambda{h}(:));
        C = Phi'*(OM.Mass{h}*U);
        for m = 1:M
            Uf = Phi*(H(:,m).*C);
            A(gOut, 1, m) = Uf;
        end
    end
end

function ef = i_find_variant(want)
    % For 'Dirac' require the real 4-block quaternion representation (Phi real, 4*nV rows)
    % that the quaternion steering math targets; legacy complex-Phi nodes are skipped.
    ef = [];
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat'));
    [~, order] = sort([d.datenum], 'descend');     % prefer the most recent matching node
    for i = order(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant');
            if ~strcmpi(m.Variant, want); continue; end
            if strcmpi(want, 'Dirac')
                E = in_bst_eigen(rel);
                hn = find(~cellfun(@isempty, E.Phi), 1);
                if isempty(hn) || ~isreal(E.Phi{hn}) || size(E.Phi{hn},1) ~= 4*numel(E.GlobalVertices{hn})
                    continue;   % skip legacy complex-Phi Dirac
                end
            end
            ef = rel; return;
        catch
        end
    end
end
