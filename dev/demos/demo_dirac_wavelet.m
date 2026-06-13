function demo_dirac_wavelet(SubjectName)
% DEMO_DIRAC_WAVELET: First experiment for the Dirac quaternion vector filterbank.
%
% Builds a Dirac eigenmode VECTOR wavelet from a vertex delta + launch direction +
% heat scale kernel, and demonstrates the CHIRALITY axis: the right-multiplication
% projectors P_+/-(n) split the field into two opposite helicities (circularly
% polarized vector wavelets rotating about n). See:
%   docs/superpowers/specs/2026-06-13-dirac-quaternion-filterbank-exploration.md
%   toolbox/forward/bst_dirac_filter.m
%
% Decisive checks (printed):
%   - no projector  -> real field, polarization helicity H = 0 (linear)
%   - P_+ and P_-   -> equal-magnitude, OPPOSITE helicity vectors (cos = -1)
%   - the helicity vector points along the chirality axis n
%
% USAGE:  demo_dirac_wavelet            % uses 'Subject01'
%         demo_dirac_wavelet('Subject01')
%
% Authors: Diellor Basha, 2026

    if (nargin < 1) || isempty(SubjectName), SubjectName = 'Subject01'; end

    % --- cortex + a lateral gyral seed vertex (global index) ---
    sSubject   = bst_get('Subject', SubjectName);
    CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
    Cx = in_tess_bst(CortexFile);
    V  = Cx.Vertices;  N = Cx.VertNormals;
    score = -V(:,1) - 2*abs(V(:,3));               % lateral (x<0), mid-height
    [~, iV] = max(score);
    Nv = N(iV,:) / norm(N(iV,:));
    a0 = [0 0 1];  dhat = a0 - (a0*Nv')*Nv;  dhat = dhat / norm(dhat);   % tangent seed
    fprintf('Seed vertex %d at [%.3f %.3f %.3f]; launch dir [%.2f %.2f %.2f]\n', iV, V(iV,:), dhat);

    % --- filter spec: heat kernel; chirality axis = ambient i (x) ---
    nModes = 400;  Tau = 0.5;
    EigTmp = tess_eigen(CortexFile, 'Dirac', 'K', nModes, 'Tau', Tau);
    h = 1;  lam = EigTmp.Lambda{h}(:);
    tHeat  = log(2) / lam(min(120, numel(lam)));   % half-attenuate ~mode 120
    Filter = struct('Type','heat','Params',struct('t',tHeat));
    nAxis  = [1 0 0];

    % --- three variants ---
    variants = {'none', 'P+', 'P-'};  signs = [0 +1 -1];
    Out = cell(1,3);  H = zeros(3,3);
    for k = 1:3
        Chir   = struct('Axis', nAxis, 'Sign', signs(k));
        Out{k} = bst_dirac_filter(CortexFile, iV, dhat, Filter, Chir, 'nModes', nModes, 'Tau', Tau);
        % polarization-helicity vector H = sum (Re V x Im V) over vertices
        H(k,:) = sum(cross(real(Out{k}.V), imag(Out{k}.V), 2), 1);
    end
    fprintf('\nvariant   polarization-helicity vector H = sum(ReV x ImV)         |H|\n');
    for k = 1:3
        fprintf('  %-5s  [%+.3e %+.3e %+.3e]   %.3e\n', variants{k}, H(k,:), norm(H(k,:)));
    end
    if norm(H(2,:)) > 0
        fprintf('P+ vs P- helicity:  cos = %.4f (expect -1)\n', (H(2,:)*H(3,:)')/(norm(H(2,:))*norm(H(3,:))));
        fprintf('no-projector |H|/|H_P+| = %.2e (expect ~0)\n', norm(H(1,:))/norm(H(2,:)));
        fprintf('helicity axis alignment (H_P+ . n)/|H_P+| = %.3f (expect ~+1)\n', (H(2,:)*nAxis(:))/norm(H(2,:)));
    end

    % --- figure: (a) wavelet arrows, (b) P+ helicity, (c) P- helicity ---
    near = find(sqrt(sum((V - V(iV,:)).^2, 2)) < 0.024);
    fpatch = Cx.Faces(all(ismember(Cx.Faces, near), 2), :);
    n=256; b=[.23 .30 .75]; w=[1 1 1]; r=[.75 .15 .15]; hh=floor(n/2);
    bwr=[linspace(b(1),w(1),hh)' linspace(b(2),w(2),hh)' linspace(b(3),w(3),hh)'; ...
         linspace(w(1),r(1),n-hh)' linspace(w(2),r(2),n-hh)' linspace(w(3),r(3),n-hh)'];

    hf = figure('Color','w','Position',[100 100 1500 540]);
    % (a) directional vector wavelet
    Va = real(Out{1}.V);  mga = sqrt(sum(Va.^2,2));
    subplot(1,3,1); hold on;
    trisurf(fpatch, V(:,1),V(:,2),V(:,3),'FaceColor',[.86 .86 .86],'EdgeColor','none','FaceAlpha',.7);
    sel = near(mga(near) > 0.2*max(mga(near)));  sc = 0.010/max(mga(sel));
    quiver3(V(sel,1),V(sel,2),V(sel,3), Va(sel,1)*sc,Va(sel,2)*sc,Va(sel,3)*sc, 0, ...
        'Color',[0 0 .7],'LineWidth',1,'MaxHeadSize',.5);
    plot3(V(iV,1),V(iV,2),V(iV,3),'r.','MarkerSize',24);
    axis equal off; view(Nv); camlight headlight; lighting gouraud;
    title('(a) vector wavelet (no chirality)','FontSize',12);
    % (b,c) helicity density
    cmax = max(abs([Out{2}.Helicity(near); Out{3}.Helicity(near)]));
    for j = 1:2
        subplot(1,3,1+j); hold on;
        trisurf(fpatch, V(:,1),V(:,2),V(:,3), Out{1+j}.Helicity, 'EdgeColor','none','FaceColor','interp','FaceLighting','none');
        plot3(V(iV,1),V(iV,2),V(iV,3),'k.','MarkerSize',24);
        axis equal off; view(Nv); caxis([-cmax cmax]); colormap(bwr);
        title(sprintf('(%c) %s helicity', 'b'+j-1, variants{1+j}),'FontSize',12);
    end
    cb = colorbar('Position',[0.93 0.3 0.012 0.4]);  cb.Label.String = '(Re V \times Im V)\cdot n';
    outpng = bst_fullfile(bst_fileparts(mfilename('fullpath')), '..', 'benchmarks', 'demo_dirac_wavelet_chirality.png');
    exportgraphics(hf, outpng, 'Resolution', 140);
    fprintf('\nWrote figure: %s\n', outpng);
end
