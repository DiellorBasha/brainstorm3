function test_covariant_operator()
% TEST_COVARIANT_OPERATOR: the 'Covariant' flat-covariant operator node passes the Hodge gates.
%
% Gate 1 (exact round-trip), Gate 3 (no curvature artifact: constant field -> 0),
% Gate 5 (complete 3-DOF). Gates 2 (matches validated div/curl) and 4 (sign) are
% exercised at the bst_helmholtz layer (test_helmholtz_covariant).
%
% USAGE:  test_covariant_operator   % Brainstorm running, TutorialAuditory loaded
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;
    sSubject = bst_get('Subject', 'Subject01');
    if isempty(sSubject) || isempty(sSubject.iCortex)
        fprintf('SKIPPED (no Subject01 cortex)\n');  fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    Cov = tess_operators(SurfaceFile, 'Covariant', 'NoSave', true);   % compute fresh, don't pollute the DB

    hh = 1;  C = Cov.Covariant{hh};
    nFh = size(C.Faces,1);  nVh = size(C.ScalarGrad,2);
    Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
    Nf = C.FaceNormal;  Af = C.FaceArea;  K = (Cov.Operator{hh}+Cov.Operator{hh}')/2;

    % smooth per-vertex scalars (spatial coordinate functions) for the in-span synthetic fields
    P = C.VertPos;  a = P(:,1);  b = P(:,3);

    % ---- GATE 1: exact round-trip on an in-span field (grad a + n x grad b) ----
    Jf = i_inspan(Gx,Gy,Gz,Nf,a,b);
    [~,~,res1] = i_hodge(Gx,Gy,Gz,Nf,Af,K,Jf);
    g1 = res1 <= 1e-10;
    fprintf('GATE1 round-trip residual = %.2e => %s\n', res1, PF{g1+1});  pass = pass && g1;

    % ---- GATE 3: constant ambient field -> machine-zero strong div/curl (even on folds) ----
    one = ones(nVh,1);  z = zeros(nVh,1);
    divC  = Gx*one + Gy*z + Gz*z;                              % strong surface divergence of [1,0,0]
    curlC = sum(cross(Nf,[Gx*z - Gy*z, Gz*one - Gx*z, Gx*z - Gz*one],2).*Nf,2); %#ok<NASGU>
    cv    = [Gy*z-Gz*z, Gz*one-Gx*z, Gx*z-Gy*one];            % strong curl vector of [1,0,0]
    omC   = sum(cv.*Nf,2);
    g3 = (max(abs(divC)) <= 1e-10) && (max(abs(omC)) <= 1e-10);
    fprintf('GATE3 constant-field max|div|=%.2e max|vort|=%.2e => %s\n', max(abs(divC)), max(abs(omC)), PF{g3+1});
    pass = pass && g3;

    % ---- GATE 5: complete 3-DOF (in-span tangential + a normal part) ----
    Jf5 = i_inspan(Gx,Gy,Gz,Nf,a,b) + 0.5*Nf;
    [~,recon5,~] = i_hodge(Gx,Gy,Gz,Nf,Af,K,Jf5);
    Jn = sum(Jf5.*Nf,2);  resid5 = Jf5 - recon5 - Jn.*Nf;
    r5 = sqrt(sum(sum(resid5.^2,2).*Af)) / sqrt(sum(sum(Jf5.^2,2).*Af));
    g5 = r5 <= 1e-10;
    fprintf('GATE5 complete 3-DOF residual = %.2e => %s\n', r5, PF{g5+1});  pass = pass && g5;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

function J = i_inspan(Gx,Gy,Gz,Nf,a,b)
    ga = [Gx*a, Gy*a, Gz*a];  gb = [Gx*b, Gy*b, Gz*b];
    J  = ga + cross(Nf, gb, 2);
end

function [phi, recon, resfrac] = i_hodge(Gx,Gy,Gz,Nf,Af,K,Jf)
    nF = size(Nf,1);  W = spdiags(Af,0,nF,nF);
    nx=Nf(:,1); ny=Nf(:,2); nz=Nf(:,3);
    Sx = spdiags(ny,0,nF,nF)*Gz - spdiags(nz,0,nF,nF)*Gy;     % (n x grad)_x
    Sy = spdiags(nz,0,nF,nF)*Gx - spdiags(nx,0,nF,nF)*Gz;
    Sz = spdiags(nx,0,nF,nF)*Gy - spdiags(ny,0,nF,nF)*Gx;
    divw  = Gx'*W*Jf(:,1) + Gy'*W*Jf(:,2) + Gz'*W*Jf(:,3);    % weak divergence
    vortw = Sx'*W*Jf(:,1) + Sy'*W*Jf(:,2) + Sz'*W*Jf(:,3);    % weak vorticity
    nV = size(K,1);  free = 2:nV;
    phi = zeros(nV,1);  psi = zeros(nV,1);
    phi(free) = K(free,free)\divw(free);
    psi(free) = K(free,free)\vortw(free);
    Virr = [Gx*phi, Gy*phi, Gz*phi];
    Vsol = cross(Nf, [Gx*psi, Gy*psi, Gz*psi], 2);
    recon = Virr + Vsol;
    resfrac = sqrt(sum(sum((Jf-recon).^2,2).*Af)) / sqrt(sum(sum(Jf.^2,2).*Af));
end
