function R = explore_spectral_pvc(SubjectName, nModes, tBlur)
% EXPLORE_SPECTRAL_PVC: exploration of spectral (Laplace-Beltrami) surface deconvolution as
% the TANGENTIAL component of partial volume correction.
%
% Idea: on the cortical sheet the PSF blur is approximated by an LBO heat low-pass exp(-t.lam)
% in the eigenbasis; PVC is the inverse. Naive inverse exp(+t.lam) amplifies noise, so we use a
% REGULARIZED (Wiener) eigenfilter  h(lam) = G/(G^2 + alpha),  G = exp(-t.lam). This sharpens
% TANGENTIALLY (within the ribbon) without touching the radial/tissue direction - complementary
% to GTM/MG. Here we measure recovery of a known cortical field vs noise, and how alpha trades
% recovery against noise amplification.
%
% Exploration only (no DB writes). One hemisphere, ico5 surface (so resolvable detail ~vertex
% spacing; a finer surface is needed for sub-mm tangential detail - noted in the takeaway).
%
% USAGE:  R = explore_spectral_pvc('sub-MTL0002', 400, 3e-5)
%
% Author: Diellor Basha, 2026

    if (nargin<1)||isempty(SubjectName), SubjectName='sub-MTL0002'; end
    if (nargin<2)||isempty(nModes), nModes=400; end
    if (nargin<3)||isempty(tBlur), tBlur=3e-5; end   % heat time of the forward tangential blur
    here = bst_fileparts(mfilename('fullpath'));

    % ----- LBO eigenbasis (tess_eigen) + mass matrix (tess_operators), one hemisphere -----
    [sS,~]=bst_get('Subject',SubjectName);
    wf = sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    fprintf('loading %d LBO eigenmodes (tess_eigen)...\n', nModes);
    EigenMat = tess_eigen(wf, 'Laplace-Beltrami', 'nModes', nModes);
    LBO = tess_operators(wf, 'Laplace-Beltrami');                 % for the mass (B) inner product
    Phi = EigenMat.Phi{1}; lam = EigenMat.Lambda{1}(:); M = LBO.Mass{1};
    vH  = double(EigenMat.GlobalVertices{1}(:));
    sSurf = in_tess_bst(wf); V = sSurf.Vertices(vH,:);
    proj = @(f) Phi' * (M * f);        % spectral coefficients (B-orthonormal basis, as bst_eigenfilter)
    rec  = @(a) Phi * a;               % reconstruct

    % ----- known cortical field with tangential structure, band-limited to the basis -----
    rng_field = 1.0 + 0.5*(V(:,1)-min(V(:,1)))/(max(V(:,1))-min(V(:,1)));   % smooth A-P gradient
    seeds = round(linspace(1,numel(vH),12)); seeds=seeds(2:end-1);
    for s=seeds, d2=sum((V-V(s,:)).^2,2); rng_field = rng_field + 1.2*exp(-d2/(2*(0.006^2))); end  % ~6mm focal blobs
    Strue = rec(proj(rng_field));      % project to the K-mode band = the recoverable signal

    % ----- forward tangential blur (heat) + transfer function -----
    G = exp(-tBlur*lam);
    Sblur = rec(G .* proj(Strue));
    sigS = std(Strue);
    fprintf('forward blur exp(-t.lam): G(lam) range [%.2f .. %.2f]; corr(blur,true)=%.3f\n', min(G),max(G), local_corr(Sblur,Strue));

    % ----- sweep noise x regularization -----
    noises = [0 0.05 0.10 0.20];           % noise SD as fraction of signal SD
    alphas = logspace(-4, 0, 25);
    R.recov = nan(numel(noises), numel(alphas));   % corr(deconv,true)
    R.base  = nan(numel(noises),1);                % corr(noisy-blur, true) no deconv
    for in=1:numel(noises)
        noise = noises(in)*sigS * local_pseudonoise(numel(vH), in);
        Sobs = Sblur + noise;
        R.base(in) = local_corr(Sobs, Strue);
        a_obs = proj(Sobs);
        for ia=1:numel(alphas)
            h = G ./ (G.^2 + alphas(ia));
            Sdec = rec(h .* a_obs);
            R.recov(in,ia) = local_corr(Sdec, Strue);
        end
    end
    R.alphas=alphas; R.noises=noises; R.lam=lam; R.G=G;

    % ----- report best alpha per noise -----
    fprintf('\n noise   no-deconv   best-deconv (alpha*)\n');
    for in=1:numel(noises)
        [bestR,ib]=max(R.recov(in,:));
        fprintf('  %4.0f%%   %7.3f     %7.3f  (a=%.1e)\n', 100*noises(in), R.base(in), bestR, alphas(ib));
    end

    % ----- figure -----
    f=figure('Visible','off','Position',[60 60 1150 430]);
    subplot(1,2,1); semilogx(R.alphas, R.recov','-','LineWidth',1.4); hold on;
    set(gca,'ColorOrderIndex',1);
    for in=1:numel(noises), yline(R.base(in),'--'); end
    xlabel('regularization \alpha'); ylabel('corr(deconvolved, true)'); ylim([0 1]); grid on;
    legend(arrayfun(@(x)sprintf('noise %.0f%%',100*x),noises,'uni',0),'Location','southwest');
    title('Spectral deconvolution recovery vs \alpha (dashed = no deconv)');
    subplot(1,2,2); plot(lam, G,'b-','LineWidth',1.4); hold on;
    for a=[1e-3 1e-2 1e-1], plot(lam, G./(G.^2+a),'-'); end
    xlabel('LBO eigenvalue \lambda'); ylabel('filter response'); grid on;
    legend({'G=exp(-t\lambda) (blur)','h, \alpha=1e-3','h, \alpha=1e-2','h, \alpha=1e-1'},'Location','northeast');
    title('Forward blur G(\lambda) and regularized inverse h(\lambda)');
    print(f, fullfile(here,'explore_spectral_pvc.png'),'-dpng','-r110'); close(f);
    fprintf('Figure: %s\n', fullfile(here,'explore_spectral_pvc.png'));
end

function v = local_pseudonoise(n, seed)
    idx=(1:n)'; v = sin(idx*(12.9898+seed)) * 43758.5453; v = 2*(v-floor(v))-1;   % deterministic [-1,1]
end

function c = local_corr(x, y)
    cc = corrcoef(x, y); c = cc(1,2);
end
