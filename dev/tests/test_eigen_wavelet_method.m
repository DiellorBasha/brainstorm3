function test_eigen_wavelet_method(ResultsFile, EigenFile)
% End-to-end: bst_eigen Method='wavelet' -> source scalogram TimefreqMat.

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

if nargin<1||isempty(ResultsFile); ResultsFile='@default_study/results_260615_1805.mat'; end
if nargin<2||isempty(EigenFile);   EigenFile='Subject01/eigen_260622_1920.mat'; end
OPTIONS = struct('Method','wavelet','EigenFile',EigenFile,'KernelName','itersine', ...
                 'Nf',8,'iTargetStudy','NoSave','TimeWindow',[202.0 202.05]);
[out, Messages, isError] = bst_eigen(ResultsFile, OPTIONS);
assert(~isError, 'bst_eigen wavelet error: %s', Messages);
FileMat = out{1};
assert(strcmpi(FileMat.Method,'eigenwavelet'), 'Method tag wrong.');
assert(ndims(FileMat.TF)==3 && size(FileMat.TF,3)>=2, 'TF must be [nSrc x nT x M].');
assert(numel(FileMat.Freqs)==size(FileMat.TF,3), 'Freqs length != #members.');
assert(all(isfinite(FileMat.TF(:))), 'TF has non-finite values.');
fprintf('Scalogram TF: [%d x %d x %d], family=%s\n', size(FileMat.TF,1), size(FileMat.TF,2), size(FileMat.TF,3), FileMat.Options.Family);
disp('test_eigen_wavelet_method PASSED');
end
