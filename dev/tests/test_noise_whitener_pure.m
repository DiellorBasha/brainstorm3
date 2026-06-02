function test_noise_whitener_pure
% Verify bst_noise_whitener (verbatim per-modality whitener from bst_inverse_linear_2018):
%   - single modality: iW whitens the covariance (iW*C*iW' ~ I)
%   - two modalities: per-modality whitening, zero cross-modality blocks
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% --- single modality (4 EEG channels): whitening property ---
rng_local = [0.7 0.2 0.1 0.0; 0.2 0.8 0.15 0.05; 0.1 0.15 0.9 0.1; 0.0 0.05 0.1 0.6];
C1 = rng_local * rng_local';                 % SPD
types1 = {'EEG','EEG','EEG','EEG'};
iW1 = bst_noise_whitener(C1, types1, 'reg', 0.1);
assert(isequal(size(iW1), [4 4]), 'iW must be [nCh x nCh].');
W = iW1 * C1 * iW1';
assert(all(isfinite(iW1(:))), 'iW must be finite.');
assert(max(abs(W - diag(diag(W))), [], 'all') < 0.5, 'Off-diagonals must be suppressed by whitening.');

% --- two modalities (2 MEG MAG + 2 EEG): cross-modality blocks zeroed ---
C2 = [2.0 0.3 0.4 0.1;
      0.3 1.7 0.2 0.2;
      0.4 0.2 0.9 0.25;
      0.1 0.2 0.25 1.1];
types2 = {'MEG MAG','MEG MAG','EEG','EEG'};
iW2 = bst_noise_whitener(C2, types2, 'reg', 0.1);
assert(max(abs(iW2(1:2,3:4)), [], 'all') == 0, 'Cross-modality whitener block must be zero.');
assert(max(abs(iW2(3:4,1:2)), [], 'all') == 0, 'Cross-modality whitener block must be zero.');
assert(all(isfinite(iW2(:))), 'iW must be finite.');

% --- 'diag' method: whitener is diagonal and finite ---
iWd = bst_noise_whitener(C1, types1, 'diag', 0.1);
assert(isequal(size(iWd), [4 4]), 'diag iW must be [nCh x nCh].');
assert(max(abs(iWd - diag(diag(iWd))), [], 'all') < 1e-12, 'diag method must yield a diagonal whitener.');
assert(all(isfinite(iWd(:))), 'diag iW must be finite.');

% --- 'none' method: no regularization, still finite and whitening ---
iWn = bst_noise_whitener(C1, types1, 'none', 0);
assert(isequal(size(iWn), [4 4]), 'none iW must be [nCh x nCh].');
assert(all(isfinite(iWn(:))), 'none iW must be finite.');
Wn = iWn * C1 * iWn';
assert(max(abs(Wn - eye(4)), [], 'all') < 1e-6, 'none method must fully whiten a full-rank covariance.');

disp('ALL TESTS PASSED');
end
