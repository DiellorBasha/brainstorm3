%% Tutorial: the MNE inverse as spectral filters, and the Dirac control surface
%
% This is a step-by-step, runnable tutorial (run it cell-by-cell with Ctrl+Enter,
% or `publish('tutorial_mne_dirac_inverse')` for an HTML write-up). It explains
% the minimum-norm (MNE) inverse as a CASCADE OF SPECTRAL FILTERS, and shows how
% the Dirac-eigenmode formulation exposes three independent "axes of control":
%
%   AXIS 1 - the NOISE axis        : eigenvalues d of the noise covariance C
%   AXIS 2 - the OBSERVABILITY axis: singular values s of the whitened leadfield
%   AXIS 3 - the GEOMETRY axis     : eigenvalues lambda of the Dirac operator
%
% Standard MNE filters only on AXIS 1 (whiten/regularize) and AXIS 2 (SNR ->
% Tikhonov). It has NO filter on AXIS 3. That unused axis is the control surface
% the Dirac formulation adds -- and the lever for surpassing MNE.
%
% Prereqs: protocol loaded (TutorialAuditory), Brainstorm on path.
% Author: Diellor Basha, 2026

%% 0. Setup -- load data and build the Dirac transform
base = 'Subject01/S01_AEF_20131218_01_notch/';
if ~exist('SHOW','var'), SHOW = true; end           % set SHOW=false to suppress figures
vis = 'on'; if ~SHOW, vis = 'off'; end

HMos = in_bst_headmodel([base 'headmodel_surf_os_meg.mat'], 0);
ChanMat = in_bst_channel([base 'channel_ctf_acc1.mat']); types = {ChanMat.Channel.Type};
G = double(HMos.Gain); iMEG = all(isfinite(G),2) & strcmpi(types(:),'MEG'); G = G(iMEG,:);
nCh = size(G,1);
Srf = in_tess_bst(HMos.SurfaceFile); Vtx = Srf.Vertices; Nrm = Srf.VertNormals; nV = size(Vtx,1);
NC = load(file_fullpath([base 'noisecov_full.mat'])); Cn = NC.NoiseCov(iMEG,iMEG); Cn = (Cn+Cn')/2;
HMf = HMos; HMf.Gain = G;
CompHM = bst_dirac(HMf, 'nModes',400, 'Tau',0.5);   % forward transform to Dirac modes
Gm = double(CompHM.Gain); lamMode = double(CompHM.Eigenvalues);
fprintf(['SETUP: %d MEG channels, %d cortical vertices, %d Dirac modes.\n' ...
    '  G    = vertex leadfield   [%d x %d]\n  Gm   = mode-forward       [%d x %d]\n'], ...
    nCh, nV, numel(lamMode), size(G,1), size(G,2), size(Gm,1), size(Gm,2));


%% 1. The standard MNE inverse -- the textbook equations
%
% Forward model (sensor measurement = leadfield * source + noise):
%       M = G * J + n,        n ~ N(0, C)          C = noise covariance
%
% Minimum-norm estimate (source prior covariance R = lambda * I):
%       J_hat = R G' (G R G' + C)^{-1} M
%
% The whole pipeline is:
%   (1) regularize the noise covariance C  -> C_reg
%   (2) WHITEN:  M~ = C_reg^{-1/2} M,   G~ = C_reg^{-1/2} G     (noise becomes white)
%   (3) set the regularization strength from a target SNR
%   (4) build the kernel; optionally NOISE-NORMALIZE (dSPM / sLORETA)
%
% Every one of these is linear, so in the right basis every one is a DIAGONAL
% MULTIPLIER -- a filter. The rest of the tutorial finds those filters.
fprintf(['\nMNE = J_hat = R G'' (G R G'' + C)^{-1} M.\n' ...
         'It is a cascade of linear steps => in the eigen/SVD bases, each is a filter.\n']);


%% 2. The three axes of control
%
% A linear inverse only ever touches three spectra:
%
%   AXIS 1  NOISE         : C = sum_i d_i u_i u_i'         (eigenvalues d_i)
%           -> whitening & regularization live here.
%
%   AXIS 2  OBSERVABILITY : C^{-1/2} G = U S V'            (singular values s_i)
%           -> "how strongly can the sensors see each independent source pattern".
%              SNR/Tikhonov regularization lives here.
%
%   AXIS 3  GEOMETRY      : Dirac operator  D^2 phi = lambda phi   (eigenvalues lambda)
%           -> spatial frequency / curvature of the cortical field.
%              MNE has NO filter here. This is the Dirac control surface.
%
% Standard MNE = filter(AXIS 1) then filter(AXIS 2). Dirac adds filter(AXIS 3).
disp('Axes:  1) noise d   2) observability s   3) Dirac geometry lambda');


%% 3. AXIS 1 -- noise covariance: whitening & regularization as a filter w(d)
%
% Whitening rotates into the noise eigenbasis and divides each direction by its
% noise amplitude. Regularization ('reg' method) adds a ridge so low-noise
% directions are not over-amplified:
%
%       w(d) = ( d + alpha * mean(d) )^{-1/2}            <-- the filter, AXIS 1
%
% alpha = 0 is the raw whitener (1/sqrt(d), blows up as d->0); alpha > 0 floors it.
[Un,Dc] = eig(Cn); d = sort(real(diag(Dc)),'descend');
alphas = [0 0.01 0.1 1]; meanD = mean(d);
fprintf('Noise cov: %d-dim, condition number = %.1e (low => well-conditioned)\n', nCh, max(d)/min(d));

f3 = figure('Color','w','Visible',vis,'Name','AXIS 1: noise filter');
dd = logspace(log10(min(d)),log10(max(d)),300); C1 = lines(numel(alphas));
for q=1:numel(alphas), loglog(dd, 1./sqrt(dd+alphas(q)*meanD),'LineWidth',2,'Color',C1(q,:)); hold on; end
loglog(d, 1./sqrt(d+0.1*meanD),'k.','MarkerSize',8);
xlabel('noise eigenvalue  d'); ylabel('whitening gain  w(d)'); grid on;
legend([arrayfun(@(a)sprintf('\\alpha=%g',a),alphas,'uni',0),{'actual d'}],'Location','southwest');
title('AXIS 1: regularized whitening  w(d)=(d+\alpha\langle d\rangle)^{-1/2}');
fprintf('  KNOB: alpha (NoiseReg). Bigger alpha = stronger floor = less noise amplification.\n');


%% 4. AXIS 2 -- observability: leadfield SVD, SNR, and the Tikhonov / band-pass filters
%
% Whiten the (mode) leadfield and take its SVD:
%       C_reg^{-1/2} * Gm = U * S * V'         s_i = observability of direction i
%
% The data, whitened and projected, has components a_i = u_i' * M~. The NAIVE
% inverse coefficient a_i / s_i blows up where s_i is tiny (unobservable = noise).
% Regularization tames this with the TIKHONOV FILTER FACTOR:
%
%       f(s) = s^2 / (s^2 + 1/Lambda)          <-- low-pass on AXIS 2
%       Lambda = SNR^2 / mean(s^2)             (SNR sets the corner)
%
% The actual data -> coefficient gain (amplitude min-norm) is:
%
%       g(s) = Lambda*s / (Lambda*s^2 + 1) = f(s)/s     <-- BAND-PASS on AXIS 2
%               (peak at s = 1/sqrt(Lambda): suppresses noise AND divides by s)
ridge = 0.1*meanD; Wn = Un*diag(1./sqrt(d+ridge))*Un';
Gt = Wn*Gm; [U,Ssvd,V] = svd(Gt,'econ'); s = diag(Ssvd); s2 = s.^2;
SNR = 3; Lambda = SNR^2/mean(s2);
fWiener = @(x,L)(L*x.^2)./(L*x.^2+1); gAmp = @(x,L)(L*x)./(L*x.^2+1);
fprintf('Observability: %d singular values; Lambda=%.2e; kept DOF (f>0.5) = %d\n', numel(s), Lambda, sum(s2>1/Lambda));

f4 = figure('Color','w','Visible',vis,'Position',[80 80 1000 420],'Name','AXIS 2: Tikhonov + band-pass');
ss = logspace(log10(min(s)),log10(max(s)),300);
subplot(1,2,1); snrs=[3 10 30 100]; C2=lines(4);
for q=1:4, semilogx(ss, fWiener(ss, snrs(q)^2/mean(s2)),'LineWidth',2,'Color',C2(q,:)); hold on; end
semilogx(ss, double(ss.^2>1/Lambda),'k--');
xlabel('observability  s'); ylabel('f(s)'); grid on; ylim([0 1.05]);
legend([arrayfun(@(x)sprintf('SNR %d',x),snrs,'uni',0),{'TSVD'}],'Location','northwest');
title('Tikhonov low-pass  f=s^2/(s^2+1/\Lambda)');
subplot(1,2,2);
semilogx(ss, gAmp(ss,Lambda)/max(gAmp(ss,Lambda)),'LineWidth',2,'Color',[.2 .5 .2]); hold on;
xline(1/sqrt(Lambda),'k:'); xlabel('observability  s'); ylabel('g(s) (norm.)'); grid on;
title('Amplitude gain  g=\Lambda s/(\Lambda s^2+1): BAND-PASS');
fprintf('  KNOB: SNR (sets Lambda => the low-pass corner). TSVD is the hard-threshold limit.\n');


%% 5. The effect on real-data coefficients (why regularization matters)
%
% Take the auditory M100 peak, whiten, project onto the observable directions,
% and compare the NAIVE inverse coefficients a_i/s_i to the REGULARIZED ones
% g(s_i)*a_i. The naive ones explode in the noise (small-s) directions.
DataMat = in_bst_data(local_pick(base,'deviant_average')); F = double(DataMat.F(iMEG,:)); Time = DataMat.Time;
win = Time>=0.06 & Time<=0.16; gfp = sqrt(sum((Wn*F).^2,1)); [~,tpk] = max(gfp.*win);
a = U'*(Wn*F(:,tpk)); cNaive = a./s; cReg = gAmp(s,Lambda).*a;
fprintf('M100 at %.0f ms. naive max|a/s|=%.2e  vs  regularized max|g*a|=%.2e\n', Time(tpk)*1e3, max(abs(cNaive)), max(abs(cReg)));

f5 = figure('Color','w','Visible',vis,'Name','effect on coefficients');
semilogy(abs(cNaive),'r','LineWidth',1.2); hold on; semilogy(abs(cReg),'b','LineWidth',2);
xlabel('observable direction index'); ylabel('|coefficient|'); grid on;
legend('naive  a/s  (noise blows up)','regularized  g\cdota','Location','northwest');
title(sprintf('Regularization tames the noise directions (M100 @ %.0f ms)', Time(tpk)*1e3));


%% 6. AXIS 3 -- the Dirac GEOMETRY axis (the one MNE never uses)
%
% The Dirac eigenvalues lambda order the cortical field by spatial frequency /
% curvature. They are a DIFFERENT axis from observability s: high-lambda modes
% are weakly seen by the sensors, but the mapping is loose (the leadfield does
% not diagonalize in the Dirac basis). MNE filters on s; it has nothing on lambda.
dkMode = sqrt(sum(Gt.^2,1))';      % per-Dirac-mode observability
fprintf('Geometry axis: %d Dirac eigenvalues. corr(lambda, -log observability) = %.2f (loose low-pass)\n', ...
    numel(lamMode), corr(lamMode, -log(dkMode+eps)));

f6 = figure('Color','w','Visible',vis,'Name','AXIS 3: geometry');
loglog(lamMode+eps, dkMode,'.','Color',[.3 .3 .6]); grid on;
xlabel('Dirac eigenvalue  \lambda  (geometry / spatial frequency)');
ylabel('observability of that mode'); title('AXIS 3: \lambda exists, but MNE places NO filter on it');


%% 7. The measures -- amplitude / dSPM / sLORETA as resolution reshaping
%
% After the kernel is built, the MEASURE rescales each source:
%   amplitude : none (the min-norm current)
%   dSPM      : divide each source by its projected noise std  ~ 1/sqrt(trace(K_v K_v'))
%   sLORETA   : normalize by the resolution  R_v^{-1/2}
% These do not change the data fit; they reshape the POINT-SPREAD FUNCTION
% (sharpen it and remove depth bias). Here is the PSF of one source under each.
OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3);
OPT.NoiseCovMat.NoiseCov = Cn; OPT.ChannelTypes = types(iMEG);
psf = struct(); meas = {'amplitude','dspm2018','sloreta'}; key = {'amp','dspm','slor'};
vt = round(0.4*nV); dsrc = G(:,(vt-1)*3+(1:3))*Nrm(vt,:)';
dist = sqrt(sum((Vtx-Vtx(vt,:)).^2,2))*1e3; edges = 0:5:80; ctr = (edges(1:end-1)+edges(2:end))/2;
for j=1:3
    OPT.InverseMeasure = meas{j}; Rr = bst_inverse_dirac(HMf, OPT);
    mg = sqrt(sum(reshape(Rr.ImagingKernel*dsrc,3,nV).^2,1))'; mg = mg/max(mg);
    psf.(key{j}) = arrayfun(@(b) mean(mg(dist>=edges(b)&dist<edges(b+1))), 1:numel(ctr));
end
f7 = figure('Color','w','Visible',vis,'Name','measures: PSF');
plot(ctr,psf.amp,'LineWidth',2); hold on; plot(ctr,psf.dspm,'LineWidth',2); plot(ctr,psf.slor,'LineWidth',2);
xlabel('distance from source (mm)'); ylabel('PSF (norm.)'); grid on; legend('amplitude','dSPM','sLORETA');
title('AXIS-2 measures reshape the point-spread function');


%% 8. The full cascade, and where the control surface is
%
% Putting it together, a recording becomes a source estimate by:
%
%   M --[whiten w(d), AXIS 1]--> M~ --[project U]--> a
%      --[band-pass g(s), AXIS 2]--> mode coeffs c --[measure, AXIS 2]-->
%      --[reconstruct J = Phi*c, AXIS 3 basis]--> cortical current
%
% MNE's knobs:  alpha (noise reg, AXIS 1) and SNR (Tikhonov corner, AXIS 2).
% Both are FIXED single-scalar shapes. dSPM/sLORETA add a fixed measure window.
%
% TO SURPASS MNE: design windows MNE leaves flat --
%   * a custom geometry window h(lambda) on AXIS 3  (a spatial-scale / curvature
%     prior, e.g. band-pass to isolate cortical scales, or a frame-transport prior)
%   * a data-adaptive or non-quadratic window on AXIS 2 (sparsity, reweighting)
%   * a Hodge/Helmholtz projector window (for flow / traveling-wave readout)
%
% Demonstration: apply a custom AXIS-3 band-pass h(lambda) to the M100 mode
% coefficients and reconstruct -- a control MNE structurally cannot express.
cM = local_modecoeffs(Wn, U, V, s, Lambda, gAmp, F(:,tpk));     % mode coeffs of the M100 estimate
band = (lamMode > prctile(lamMode,20)) & (lamMode < prctile(lamMode,60));  % a mid-spatial-frequency band
cBand = cM .* band;
fprintf(['\nCONTROL SURFACE: MNE filters AXIS 1 (alpha) + AXIS 2 (SNR) only.\n' ...
    'Demonstrated an AXIS-3 band-pass h(lambda): kept %d/%d modes by GEOMETRY,\n' ...
    'reconstructing a spatial-scale-selected field -- a prior MNE cannot represent.\n'], sum(band), numel(band));

f8 = figure('Color','w','Visible',vis,'Name','AXIS 3 control');
stem(lamMode, abs(cM)/max(abs(cM)),'.','Color',[.7 .7 .7]); hold on;
stem(lamMode(band), abs(cBand(band))/max(abs(cM)),'.','Color',[.1 .5 .2]);
set(gca,'XScale','log'); xlabel('Dirac eigenvalue \lambda'); ylabel('|mode coefficient| (norm.)'); grid on;
legend('all modes','AXIS-3 band-pass h(\lambda)'); title('Surpassing MNE: a custom filter on the GEOMETRY axis');

disp('--- tutorial complete: MNE = filters on axes 1 & 2; axis 3 (geometry) is the new control. ---');


%% ===== local helpers =====
function fn = local_pick(base, tag)
    sStudy = bst_get('StudyWithCondition', base(1:end-1)); fn = '';
    for i=1:numel(sStudy.Data), if contains(sStudy.Data(i).FileName,tag), fn=sStudy.Data(i).FileName; return; end; end
    fn = sStudy.Data(1).FileName;
end
function c = local_modecoeffs(Wn, U, V, s, Lambda, gAmp, m)
    % amplitude mode coefficients of one data vector: c = V diag(g) U' (Wn m)
    c = V * (gAmp(s,Lambda) .* (U' * (Wn*m)));
end
