function test_bst_divergence_ambient
% TEST_BST_DIVERGENCE_AMBIENT: Parity test for the ambient bst_divergence branch.
% Verifies that the flat-covariant surface divergence matches the helmholtz_baseline.mat oracle.
%
% USAGE:
%   test_bst_divergence_ambient
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@

    d = load(fullfile(fileparts(mfilename('fullpath')), 'baselines', 'helmholtz_baseline.mat'));
    B = d.B;  Cov = tess_operators(B.Surf, 'Covariant');
    div = bst_divergence(B.J, [], 'Ambient', [], Cov);
    assert(isequal(size(div), size(B.Div)), 'div shape mismatch');
    rel = norm(div - B.Div) / max(norm(B.Div), eps);
    assert(rel < 1e-10, 'ambient divergence differs from baseline (rel=%.2e)', rel);
    fprintf('PASS test_bst_divergence_ambient (rel=%.2e)\n', rel);
end
