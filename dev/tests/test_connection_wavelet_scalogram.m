function test_connection_wavelet_scalogram()
% TEST_CONNECTION_WAVELET_SCALOGRAM: the bst_eigen Method='wavelet' FILE path produces a
% 3-D-tangent scalogram for the Connection Laplacian from a real 3-vector source map, at
% parity with the scalar LBO and Dirac paths.
%
% Validates the tangent<->complex frame bridge (GetConnectionFrame/TangentToComplex/
% ComplexToTangent in bst_eigen):
%   - tangent -> complex -> tangent is a machine-exact isometry on the tangent part;
%   - the full pipeline (encode -> Analysis -> Synthesis -> decode) reconstructs the
%     band-limited tangent field exactly (tight frame);
%   - bst_eigen(F3, Method='wavelet') on a connection node returns a real 3-vector
%     scalogram [3nV x nT x M], nComponents=3, Method='eigenwavelet'.
%
% Needs nxr-compute + a Connection Laplacian eigen node under Subject01. SKIPs if absent.
%
% Authors: Diellor Basha, 2026

    nFail = 0;  chk = @(n,c) i_chk(n,c);
    efC = i_find_variant('Connection Laplacian');
    if isempty(efC)
        fprintf('SKIP test_connection_wavelet_scalogram: no Connection Laplacian eigen node.\n'); return;
    end
    fprintf('Connection node: %s\n', efC);

    EC = in_bst_eigen(efC);  OC = in_bst_operator(EC.OperatorFile);
    nVg = max(cellfun(@(v) max([v(:);0]), EC.GlobalVertices));
    nT = 2;

    % global frame (same as bst_eigen's GetConnectionFrame)
    E1 = zeros(nVg,3); E2 = zeros(nVg,3);
    for h = 1:numel(EC.GlobalVertices)
        gv = EC.GlobalVertices{h}(:); if isempty(gv); continue; end
        Fr = bst_operator_frame(OC, h); E1(gv,:) = Fr.e1; E2(gv,:) = Fr.e2;
    end

    % in-subspace complex field -> decode to a PURE-TANGENT 3-vector field F3
    Z0 = zeros(nVg, nT);
    for h = 1:numel(EC.Phi)
        P = EC.Phi{h}; if isempty(P); continue; end
        Z0(EC.GlobalVertices{h},:) = P*(randn(size(P,2),nT) + 1i*randn(size(P,2),nT));
    end
    F3 = i_c2t(Z0, E1, E2);

    % ---- (1) tangent encode/decode is a machine-exact isometry ----
    Z1 = i_t2c(F3, E1, E2);
    nFail = nFail + chk('tangent <-> complex isometry (exact)', norm(Z1(:)-Z0(:))/max(norm(Z0(:)),eps) < 1e-12);

    % ---- (2) full pipeline reconstructs the band-limited tangent field (tight frame) ----
    fr = bst_eigenwavelet('Design', 'itersine', 6, i_grange(EC));
    Fc = i_t2c(F3, E1, E2);
    Wc = bst_eigenwavelet('Analysis', Fc, EC, OC, fr);
    Zrec = bst_eigenwavelet('Synthesis', Wc, EC, OC, fr);
    % band-limited complex projection (in-subspace => equals Fc)
    proj = zeros(nVg, nT);
    for h = 1:numel(EC.Phi)
        P = EC.Phi{h}; if isempty(P); continue; end; gv = EC.GlobalVertices{h}(:);
        proj(gv,:) = P*(P'*(OC.Mass{h}*Fc(gv,:)));
    end
    nFail = nFail + chk('tight-frame pipeline RT == projection (exact)', norm(Zrec(:)-proj(:))/max(norm(proj(:)),eps) < 1e-9);
    F3rec = i_c2t(Zrec, E1, E2);
    nFail = nFail + chk('decoded 3-vector RT == input tangent field', norm(F3rec(:)-F3(:))/max(norm(F3(:)),eps) < 1e-9);

    % ---- (3) bst_eigen FILE path: real 3-vector -> 3-vector scalogram ----
    O = bst_eigen(); O.Method = 'wavelet'; O.EigenFile = efC; O.KernelName = 'itersine'; O.Nf = 6;
    [out, msg, err] = bst_eigen(F3, O);
    nFail = nFail + chk('bst_eigen wavelet: no error', err == 0 && isempty(msg));
    if iscell(out) && ~isempty(out)
        R = out{1}; if ischar(R); R = load(R); end
        nFail = nFail + chk('scalogram TF [3nV x nT x M], real', isequal(size(R.TF),[3*nVg, nT, numel(fr.g)]) && isreal(R.TF));
        nFail = nFail + chk('scalogram nComponents==3', R.nComponents == 3);
        nFail = nFail + chk('scalogram Method==eigenwavelet', strcmp(R.Method,'eigenwavelet'));
    else
        nFail = nFail + chk('bst_eigen returned a scalogram', false);
    end

    fprintf('\n==== test_connection_wavelet_scalogram: %d failed ====\n', nFail);
    if nFail > 0; error('test_connection_wavelet_scalogram FAILED'); end
    disp('ALL TESTS PASSED');
end

% ===== helpers (mirror bst_eigen's private bridge) =====
function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end
function Z = i_t2c(F3, E1, E2)
    nVg = size(E1,1); nT = size(F3,2); F3r = reshape(F3, 3, nVg, nT);
    Z = reshape(sum(F3r.*reshape(E1.',3,nVg),1), nVg, nT) + 1i*reshape(sum(F3r.*reshape(E2.',3,nVg),1), nVg, nT);
end
function V = i_c2t(Z, E1, E2)
    [nVg, nT, M] = size(Z); V = zeros(3*nVg, nT, M);
    for m = 1:M
        Zm = Z(:,:,m); Vk = zeros(3, nVg, nT);
        for k = 1:3; Vk(k,:,:) = real(Zm).*E1(:,k) + imag(Zm).*E2(:,k); end
        V(:,:,m) = reshape(Vk, 3*nVg, nT);
    end
end
function r = i_grange(E)
    lams = E.Lambda(~cellfun(@isempty, E.Lambda)); r = [min(cellfun(@min,lams)), max(cellfun(@max,lams))];
end
function ef = i_find_variant(want)
    ef = [];
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat')); [~,ord] = sort([d.datenum],'descend');
    for i = ord(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant');
            if strcmpi(m.Variant, want); ef = rel; return; end
        catch
        end
    end
end
