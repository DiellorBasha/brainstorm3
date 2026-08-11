function demo_cortical_flow_core(SurfaceFile)
% DEMO_CORTICAL_FLOW_CORE: GUI-free tour of the cortical-flow core functions
% (tess_massmatrix, tess_laplacian, tess_operators, tess_eigen) on the
% feature/cortical-flow-core branch.
%
% HOW TO RUN (the clean branch must NOT use your real ~/.brainstorm — your
% dev-migrated config triggers a DB-version dialog):
%     ./dev/launch_clean.sh              % opens desktop MATLAB + isolated Brainstorm
% then in that MATLAB:
%     demo_cortical_flow_core            % uses the Gate-0 sub-0002 ico5 cortex
%     demo_cortical_flow_core(file)      % or any surface .mat with a Structures atlas
%
% All figures are plain MATLAB (trisurf) — no Brainstorm viewers involved.

if nargin < 1 || isempty(SurfaceFile)
    SurfaceFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'dev', 'verify', 'phase0', 'bst_userdir_clean', '.brainstorm', 'local_db', ...
        'omega-tutorial-cortical-flow', 'anat', 'sub-0002', 'tess_cortex_pial_low.mat');
end
assert(exist(SurfaceFile, 'file') == 2, 'Surface not found: %s', SurfaceFile);

%% 1) Primitives on an analytic sphere ------------------------------------
fprintf('=== 1) Primitives: unit icosphere (2562 vertices) ===\n');
[Vs, Fs] = tess_sphere(2562);
Bs = tess_massmatrix(Vs, Fs);
As = tess_laplacian(Vs, Fs);
fprintf('  total mass  = %.6f   (sphere area 4*pi = %.6f)\n', full(sum(Bs(:))), 4*pi);
fprintf('  ||A*1||_inf = %.3g   (constants in null space)\n', max(abs(As*ones(size(As,1),1))));
fprintf('  symmetry    : A %.3g, B %.3g (Frobenius asymmetry)\n', ...
    norm(As-As','fro'), norm(Bs-Bs','fro'));

%% 2) Pencil on the real cortex -------------------------------------------
fprintf('\n=== 2) tess_operators: per-hemisphere LBO pencil ===\n');
[Op, Ms, vH] = tess_operators(SurfaceFile, 'Laplace-Beltrami');
for hh = 1:2
    fprintf('  hemi %d (%s): %d vertices, stiffness nnz %d, mass sum %.4g m^2\n', ...
        hh, char('L'+(hh-1)*('R'-'L')), size(Op{hh},1), nnz(Op{hh}), full(sum(Ms{hh}(:))));
end

%% 3) Eigenbasis: compute (or cache-hit), inspect, reuse -------------------
fprintf('\n=== 3) tess_eigen: 400 modes/hemisphere (shift-invert) ===\n');
t0 = tic; E = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', 400); t1 = toc(t0);
fprintf('  first call : %.2f s  (%s)\n', t1, ternary(t1 < 1, 'served from embedded cache', 'fresh eigensolve'));
t0 = tic; E = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', 400); t2 = toc(t0);
fprintf('  second call: %.2f s  (always cached)\n', t2);
fprintf('  solver: %s, TauRel=%g | lambda_1 = [%.2g, %.2g], lambda_400 = [%.5g, %.5g]\n', ...
    E.Solver.Method, E.Solver.TauRel, E.Lambda{1}(1), E.Lambda{2}(1), E.Lambda{1}(400), E.Lambda{2}(400));

figure('Name', 'LBO spectra', 'Color', 'w');
semilogy(E.Lambda{1}, 'LineWidth', 1.2); hold on; semilogy(E.Lambda{2}, 'LineWidth', 1.2);
grid on; xlabel('mode index'); ylabel('\lambda  [m^{-2}]');
legend({'left', 'right'}, 'Location', 'southeast'); title('LBO eigenvalue spectra per hemisphere');

% spatial wavelength calibration: 2*pi/sqrt(lambda), in mm
fprintf('  wavelength 2\\pi/sqrt(lambda): mode 10 ~ %.0f mm, mode 100 ~ %.0f mm, mode 400 ~ %.0f mm\n', ...
    2*pi/sqrt(E.Lambda{1}(10))*1000, 2*pi/sqrt(E.Lambda{1}(100))*1000, 2*pi/sqrt(E.Lambda{1}(400))*1000);

%% 4) Eigenmodes + spectral (heat-kernel) smoothing, rendered with trisurf --
fprintf('\n=== 4) Modes + heat-kernel smoothing (plain trisurf render) ===\n');
T = load(SurfaceFile);
Faces = double(T.Faces); nVtot = size(T.Vertices, 1);
figure('Name', 'Eigenmodes + heat smoothing', 'Color', 'w');
modes = [3 10 50];
for m = 1:numel(modes)
    subplot(2, numel(modes), m);
    PlotBothHemis(T.Vertices, Faces, vH, {E.Phi{1}(:,modes(m)), E.Phi{2}(:,modes(m))});
    title(sprintf('mode %d  (\\lambda=%.0f m^{-2})', modes(m), E.Lambda{1}(modes(m))));
end
% heat smoothing of a random spike map on the LEFT hemisphere: y = Phi exp(-t lambda) Phi' B x
rng(1);
xL = zeros(size(Op{1},1), 1); xL(randperm(numel(xL), 25)) = 1;   % 25 random spikes
cL = E.Phi{1}' * (Ms{1} * xL);                                    % analysis (B-inner product)
tScales = [0, 4/E.Lambda{1}(200), 4/E.Lambda{1}(30)];             % none, mid, heavy
for k = 1:3
    yL = E.Phi{1} * (exp(-tScales(k) * E.Lambda{1}) .* cL);       % synthesis
    subplot(2, 3, 3+k);
    PlotBothHemis(T.Vertices, Faces, vH, {yL, zeros(size(Op{2},1),1)});
    title(sprintf('heat smoothing t=%.3g', tScales(k)));
end
fprintf('  done — two figures should be open.\n');
end

function PlotBothHemis(Vertices, Faces, vH, valsPerHemi)
    vals = nan(size(Vertices,1), 1);
    for hh = 1:2, vals(vH{hh}) = valsPerHemi{hh}; end
    isV = false(size(Vertices,1),1); isV([vH{1}; vH{2}]) = true;
    fMask = all(isV(Faces), 2);
    trisurf(Faces(fMask,:), Vertices(:,1), Vertices(:,2), Vertices(:,3), vals, ...
        'EdgeColor', 'none');
    axis equal off vis3d; view(0, 90); camlight headlight; lighting gouraud;
    colormap(gca, turbo); colorbar('southoutside');
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
