% EXPERIMENT_SHIFT_INVERT: pin the LBO shift-invert parameters (Open Question 1).
% Arms: shift-invert at TauRel in {1e-6, 1e-4, 1e-2} (run FIRST, safe) vs
% sigma=0 'smallestabs' (comparison only, run LAST/isolated -- this arm is
% known to be fatally ill-conditioned on the real cortex pencil and can kill
% the MATLAB process outright, not just throw a catchable error). K=400.
% Metrics: residual, B-orthonormality, lambda_1 zero-mode recovery,
% run-to-run reproducibility (2 runs / arm, principal-angle subspace
% correlation), cross-arm subspace agreement, runtime.
%
% Invocation is split two ways for crash isolation (see 2026-08-08 rerun):
%   env PHASE2_HEMI = '1' | '2'   (required)
%   env PHASE2_ARM  = 'shift' | 'sigma0'   (required)
% 'shift' runs the 3 shift-invert arms + their pairwise cross-arm agreement
% and APPENDS results incrementally (open/write/close per arm) so a crash
% mid-run loses at most the in-flight arm, never prior arms. 'sigma0' runs
% only the sigma=0 comparison arm, in its own process, so a fatal death there
% cannot take the shift-invert results down with it.
resultsFile = fullfile(fileparts(mfilename('fullpath')), 'shift_experiment_results.md');
S = load(fullfile(fileparts(mfilename('fullpath')), '..', 'phase1', 'oracle_lbo_sub0002.mat'));
K = 400;
TauRelList = [1e-6, 1e-4, 1e-2];

hemiStr = getenv('PHASE2_HEMI');
armStr  = getenv('PHASE2_ARM');
assert(any(strcmp(hemiStr, {'1','2'})), 'PHASE2_HEMI must be ''1'' or ''2''');
assert(any(strcmp(armStr, {'shift','sigma0'})), 'PHASE2_ARM must be ''shift'' or ''sigma0''');
hh = str2double(hemiStr);

A = S.A{hh}; B = S.B{hh}; n = size(A,1);
sigmaScale = full(sum(diag(A))) / full(sum(diag(B)));

% Header: only the hemi==1 + arm=='shift' invocation creates the file fresh;
% every other invocation appends.
isFirstInvocation = (hh == 1) && strcmp(armStr, 'shift');
fid = fopen(resultsFile, ternary(isFirstInvocation, 'w', 'a'));
if isFirstInvocation
    fprintf(fid, '# Shift-invert experiment (K=%d, real cortex pencil)\n\n', K);
end
if strcmp(armStr, 'shift')
    fprintf(fid, '## Hemisphere %d (n=%d, sigmaScale=%.6g)\n\n', hh, n, sigmaScale);
end
fclose(fid);

if strcmp(armStr, 'shift')
    % --- shift-invert arms only, safe, run first, write durably per arm ---
    armPhi = cell(1, numel(TauRelList));
    for arm = 1:numel(TauRelList)
        label = sprintf('shift-invert TauRel=%g', TauRelList(arm));
        tau = TauRelList(arm) * sigmaScale;
        try
            t0 = tic;
            [Phi1, L1] = eigsArm(A, B, K, tau, ones(n,1));
            t1 = toc(t0);
            [Phi2, L2] = eigsArm(A, B, K, tau, ones(n,1));  %#ok<ASGLU> % identical start: determinism probe
            [Phi3, L3] = eigsArm(A, B, K, tau, rand(n,1));  %#ok<ASGLU> % different start: robustness probe
            res  = pencilResidual(A, B, Phi1, L1);
            orth = norm(Phi1' * B * Phi1 - eye(K), 'fro');
            rep12 = subspaceCorr(Phi1, Phi2, B);
            rep13 = subspaceCorr(Phi1, Phi3, B);
            line = sprintf('- %s: t=%.1fs, maxres=%.3g, orth=%.3g, lambda1=%.3g, det-rep=%.3g, rand-rep=%.3g\n', ...
                label, t1, res, orth, L1(1), 1-rep12, 1-rep13);
            armPhi{arm} = Phi1;
        catch err
            line = sprintf('- %s: FAILED (%s)\n', label, err.message);
            armPhi{arm} = [];
        end
        fid = fopen(resultsFile, 'a'); fprintf(fid, '%s', line); fclose(fid);
    end
    % --- cross-arm agreement among the shift-invert arms ---
    for a = 1:numel(armPhi)
        for b = a+1:numel(armPhi)
            if ~isempty(armPhi{a}) && ~isempty(armPhi{b})
                line = sprintf('- cross-arm subspace disagreement arms %d vs %d: %.3g\n', ...
                    a, b, 1 - subspaceCorr(armPhi{a}, armPhi{b}, B));
                fid = fopen(resultsFile, 'a'); fprintf(fid, '%s', line); fclose(fid);
            end
        end
    end
    fid = fopen(resultsFile, 'a'); fprintf(fid, '\n'); fclose(fid);
    disp('SHIFT ARMS DONE');
else
    % --- comparison arm: plain smallestabs, sigma=0 (known-bug arm) ---
    % Isolated invocation: if this dies fatally (not a catchable MATLAB
    % error -- e.g. process death during an ill-conditioned factorization),
    % it takes only this process down, not the shift-invert results already
    % written to disk.
    label = 'sigma0-smallestabs';
    try
        t0 = tic;
        [Phi1, L1] = eigsArm(A, B, K, 0, ones(n,1));
        t1 = toc(t0);
        [Phi2, L2] = eigsArm(A, B, K, 0, ones(n,1));  %#ok<ASGLU>
        [Phi3, L3] = eigsArm(A, B, K, 0, rand(n,1));  %#ok<ASGLU>
        res  = pencilResidual(A, B, Phi1, L1);
        orth = norm(Phi1' * B * Phi1 - eye(K), 'fro');
        rep12 = subspaceCorr(Phi1, Phi2, B);
        rep13 = subspaceCorr(Phi1, Phi3, B);
        line = sprintf('- %s: t=%.1fs, maxres=%.3g, orth=%.3g, lambda1=%.3g, det-rep=%.3g, rand-rep=%.3g\n', ...
            label, t1, res, orth, L1(1), 1-rep12, 1-rep13);
    catch err
        line = sprintf('- %s: FAILED (%s)\n', label, err.message);
    end
    fid = fopen(resultsFile, 'a'); fprintf(fid, '%s', line); fclose(fid);
    disp('SIGMA0 ARM DONE');
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
function [Phi, Lambda] = eigsArm(A, B, K, tau, startVec)
    % R2023b eigs accepts the legacy opts-struct form (opts.v0/opts.issym);
    % fall back to the name-value 'StartVector' form if that errors on a
    % different MATLAB version. issym is irrelevant for matrix inputs (eigs
    % detects symmetry automatically); it is set here only as a documented
    % hint in the legacy-opts path.
    try
        opts.v0 = startVec; opts.issym = true;
        [Phi, Mu] = eigs(A + tau*B, B, K, 'smallestabs', opts);
    catch
        [Phi, Mu] = eigs(A + tau*B, B, K, 'smallestabs', 'StartVector', startVec);
    end
    [Lambda, iSort] = sort(diag(Mu) - tau);
    Phi = Phi(:, iSort);
end
function r = pencilResidual(A, B, Phi, Lambda)
    R = A*Phi - B*Phi*diag(Lambda);
    N = (A + B)*Phi;
    r = max(sqrt(sum(R.^2,1)) ./ sqrt(sum(N.^2,1)));
end
function c = subspaceCorr(P1, P2, B)
    % mean squared cosine of principal angles between B-orthonormal bases
    M = P1' * B * P2;
    c = mean(svd(M).^2);
end
