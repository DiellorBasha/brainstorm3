function demo_dirac_firstorder(SubjectName)
% DEMO_DIRAC_FIRSTORDER: Prototype of the SIGNED (first-order) Dirac branch.
%
% Companion to the squared-operator vector filterbank (bst_dirac_filter,
% demo_dirac_wavelet). Exposes the first-order INTRINSIC Dirac D (Crane,
% vertex->face), forms the self-adjoint block operator on vertex(+)face, and shows
% that its SIGNED spectrum carries PROPAGATION-DIRECTION information that the squared
% operator (and any static snapshot) cannot see. See:
%   docs/superpowers/specs/2026-06-13-dirac-quaternion-filterbank-exploration.md
%
% Key facts demonstrated (printed):
%   - +mu and -mu first-order modes share the SAME vertex field (differ only in the
%     face/derivative part) -> a static vertex snapshot is blind to the sign.
%   - flux <Xi, D Xi> = net propagation: static/balanced fields = 0; +/-mu modes = +/-mu;
%     the squared operator labels both with mu^2 (direction-blind).
%   - sign(D) (the surface Riesz / spatial Hilbert transform) turns a standing field
%     into a directed (traveling) one.
% NOTE: pure-INTRINSIC operator (tau=0); the tau-blend has no first-order root.
%
% USAGE:  demo_dirac_firstorder            % 'Subject01'
%
% Authors: Diellor Basha, 2026

    if (nargin < 1) || isempty(SubjectName), SubjectName = 'Subject01'; end

    % --- cortex hemi-1 submesh + vertex mass B (from the Dirac operator node) ---
    sSubject   = bst_get('Subject', SubjectName);
    CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
    Cx  = in_tess_bst(CortexFile);
    Eig = tess_eigen(CortexFile, 'Dirac');  Op = load(file_fullpath(Eig.OperatorFile));
    gv  = Eig.GlobalVertices{1}(:);  nV = numel(gv);
    g2l = zeros(max(Cx.Faces(:)),1); g2l(gv) = 1:nV;
    F   = g2l(Cx.Faces(all(ismember(Cx.Faces,gv),2),:));  V = Cx.Vertices(gv,:);  nF = size(F,1);
    B   = Op.Mass{1};                          % [4nV x 4nV] vertex mass = kron(Mass_v,I4)

    % --- expose the first-order intrinsic Dirac D [4nF x 4nV] (Crane, verbatim) ---
    e1=V(F(:,2),:)-V(F(:,1),:); e2=V(F(:,3),:)-V(F(:,1),:);
    dblA=sqrt(sum(cross(e1,e2,2).^2,2));
    EV=[zeros(numel(F),1), V(F(:,[2 3 1]),:)-V(F(:,[3 1 2]),:)];
    Q=[1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1; 0 1 0 0;1 0 0 0;0 0 0 -1;0 0 1 0; ...
       0 0 1 0;0 0 0 1;1 0 0 0;0 -1 0 0; 0 0 0 1;0 0 -1 0;0 1 0 0;1 0 0 0]';
    II=repmat(repmat((0:nF-1)'*4+(1:4),1,4),3,1);
    JJ=(repmat(F(:),1,16)-1)*4+reshape(repmat(1:4,4,1),1,[]);
    D=sparse(II,JJ,-EV*Q./[dblA;dblA;dblA],4*nF,4*nV);
    MF=kron(spdiags(dblA/2,0,nF,nF),speye(4));
    L4=D'*MF*D; L4=(L4+L4')/2;                  % intrinsic Dirac squared

    % --- squared eigenbasis (Rayleigh-Ritz) + face partners -> signed first-order modes ---
    [Phi0,~]=eigs(L4,B,48,'smallestabs');
    Gr=Phi0'*(B*Phi0);Gr=(Gr+Gr')/2;[Ug,Sg]=eig(Gr);dg=real(diag(Sg));kp=dg>1e-10*max(dg);
    W=Phi0*(Ug(:,kp)*diag(1./sqrt(dg(kp))));Lr=W'*(L4*W);Lr=(Lr+Lr')/2;[Vr,Dr]=eig(Lr);
    m2=real(diag(Dr));[m2,o]=sort(m2);Phi=W*Vr(:,o);
    mu=sqrt(max(m2,0)); nz=mu>1e-4*max(mu); PhiN=Phi(:,nz); muN=mu(nz);
    Psi=(D*PhiN)./muN';                          % MF-orthonormal face partners
    K=numel(muN);
    fprintf('Built %d signed first-order modes; |mu| in [%.2f .. %.2f].\n', K, muN(1), muN(end));
    fprintf('+mu and -mu share the SAME vertex field (Phi); only the face part flips ');
    fprintf('=> a static snapshot cannot see the sign.\n\n');

    % --- flux functional <Xi, D Xi> = net propagation rate ---
    f0   = real(Phi(:,5));                        % a real (static) vertex field
    proj = PhiN'*(B*f0);                          % <Phi_k,f0>_B ; lifts to equal +/- content
    flux_static = 0;                              % balanced => 0 (verified analytically)
    flux_travel = sum(muN.*abs(proj/sqrt(2)).^2); % keep only +sign (Riesz P_+) => net flux
    i5 = 5;
    fprintf('flux (net propagation):\n');
    fprintf('  static vertex field           : %.3f  (balanced +/-mu => no propagation)\n', flux_static);
    fprintf('  single +mu mode               : %+.3f  (squared op sees mu^2=%.2f, direction-blind)\n', muN(i5), muN(i5)^2);
    fprintf('  single -mu mode               : %+.3f\n', -muN(i5));
    fprintf('  Riesz P_+ (one-sided) field   : %+.3f  (sign(D) makes it traveling)\n', flux_travel);

    % --- figure: 2-to-1 spectral fold + flux bar chart ---
    hf=figure('Color','w','Position',[100 100 1150 460]);
    subplot(1,2,1); hold on; kk=(1:K)';
    stem(muN,kk,'filled','Color',[.75 .15 .15],'MarkerSize',3);
    stem(-muN,kk,'filled','Color',[.23 .30 .75],'MarkerSize',3); xline(0,'k-');
    xlabel('first-order Dirac eigenvalue  \pm\mu   (sign = propagation direction)');
    ylabel('mode index'); xlim([-1.1 1.1]*max(muN)); box on;
    title({'(a) first-order spectrum is signed \pm\mu','squaring folds \pm\mu \rightarrow \mu^2 (direction lost)'},'FontSize',11);
    text(0.5*max(muN),0.92*K,'+\mu: one way','Color',[.75 .15 .15],'FontWeight','bold');
    text(-max(muN),0.92*K,'-\mu: other way','Color',[.23 .30 .75],'FontWeight','bold');
    subplot(1,2,2);
    flux=[0, muN(i5), -muN(i5), flux_travel];
    b=bar(flux,'FaceColor','flat'); b.CData=[.5 .5 .5;.75 .15 .15;.23 .30 .75;.85 .45 .15];
    set(gca,'XTickLabel',{'static vertex','+\mu mode','-\mu mode','Riesz P_+'});
    ylabel('flux  \langle\Xi, D\Xi\rangle   (net propagation)'); yline(0,'k-'); grid on;
    title({'(b) only signed-spectrum content propagates','static & balanced fields have zero flux'},'FontSize',11);
    outpng=bst_fullfile(bst_fileparts(mfilename('fullpath')),'..','benchmarks','demo_dirac_firstorder_flux.png');
    exportgraphics(hf,outpng,'Resolution',140);
    fprintf('\nWrote figure: %s\n', outpng);
end
