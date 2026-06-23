function test_eigenfilter_orchestrator(EigenFile)
% Verb dispatch + Analysis correctness vs an independent manifold computation.

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
%
% Authors: Diellor Basha, 2026

if nargin < 1 || isempty(EigenFile); EigenFile = 'Subject01/eigen_260622_1920.mat'; end
E = in_bst_eigen(EigenFile);
O = in_bst_operator(E.OperatorFile);
% scalar Laplace-Beltrami map: one row per global vertex
nSrc = max(cellfun(@(v) max([v(:);0]), E.GlobalVertices));
rng(7); F = randn(nSrc, 5);
% --- Design + Evaluate ---
g = bst_eigenfilter('Design', 'heat', struct('t',0.02));
h = bst_eigenfilter('Evaluate', g, E.Lambda{1}(:));
assert(numel(h)==numel(E.Lambda{1}), 'Evaluate gain length mismatch.');
% --- Analysis vs independent math: Phi*(h.*(Phi'*(B*F_h))) per hemisphere ---
Ffilt = bst_eigenfilter('Analysis', F, E, O, 'heat', struct('t',0.02));
Ref = zeros(size(F));
for hh=1:numel(E.Phi)
    Phi=E.Phi{hh}; if isempty(Phi); continue; end
    idx=E.GlobalVertices{hh}(:); B=O.Mass{hh}; lam=E.Lambda{hh}(:);
    hh_gain = exp(-0.02*lam);
    Ref(idx,:) = Phi*(hh_gain .* (Phi'*(B*F(idx,:))));
end
assert(norm(Ffilt(:)-Ref(:))/norm(Ref(:)) < 1e-10, 'Analysis != independent manifold filter.');
% --- RowMap ---
[sr,dr,nr,msg] = bst_eigenfilter('RowMap', F, E, 1); %#ok<ASGLU>
assert(isempty(msg) && nr>=numel(E.GlobalVertices{1}), 'RowMap failed.');
disp('test_eigenfilter_orchestrator PASSED');
end
