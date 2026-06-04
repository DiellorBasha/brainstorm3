function TessMat = out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Vertices, Faces, isInteractive)
% OUT_TESS_CONN_EIGENMODES: Save connection-Laplacian eigenmodes into a surface file.
%
% USAGE:  TessMat = out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Vertices, Faces)
%         TessMat = out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Vertices, Faces, isInteractive)
%
% SEE ALSO: in_tess_conn_eigenmodes, tess_conn_eigenmodes

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

if nargin < 5
    isInteractive = true;
end

%% ===== VALIDATE INPUTS =====
if ~isstruct(ConnEig) || ~isfield(ConnEig, 'Vectors') || ~isfield(ConnEig, 'Values')
    error('ConnEig must be a structure with Vectors and Values fields.');
end
if size(ConnEig.Vectors, 1) ~= size(Vertices, 1)
    error('ConnEig.Vectors (%d rows) must match Vertices (%d rows).', ...
        size(ConnEig.Vectors, 1), size(Vertices, 1));
end

%% ===== LOAD SURFACE FILE =====
SurfaceFileFull = file_fullpath(SurfaceFile);
if ~file_exist(SurfaceFileFull)
    error('Surface file not found: %s', SurfaceFileFull);
end
TessMat = load(SurfaceFileFull);

%% ===== STORE (complex single vectors; double sparse operators) =====
Store = struct();
Store.Vectors     = single(ConnEig.Vectors);     % complex single
Store.Values      = ConnEig.Values(:);
Store.nModes      = ConnEig.nModes;
if isfield(ConnEig, 'Component'),   Store.Component   = ConnEig.Component(:);   end
if isfield(ConnEig, 'CompRank'),    Store.CompRank    = ConnEig.CompRank(:);    end
if isfield(ConnEig, 'Order'),       Store.Order       = ConnEig.Order(:);       end
if isfield(ConnEig, 'nComponents'), Store.nComponents = ConnEig.nComponents;    end
% Sparse stays double (MATLAB has no single sparse); negligible vs Vectors.
if isfield(ConnEig, 'MassMatrix')    && ~isempty(ConnEig.MassMatrix),    Store.MassMatrix    = ConnEig.MassMatrix;    end
if isfield(ConnEig, 'ConnLaplacian') && ~isempty(ConnEig.ConnLaplacian), Store.ConnLaplacian = ConnEig.ConnLaplacian; end
Store.OperatorType   = ConnEig.OperatorType;
Store.nSym           = ConnEig.nSym;
Store.Regularization = ConnEig.Regularization;
Store.Sigma          = ConnEig.Sigma;
Store.Tolerance      = ConnEig.Tolerance;
Store.nRemoved       = ConnEig.nRemoved;
Store.ComputeTime    = ConnEig.ComputeTime;
Store.ComputeDate    = datestr(now, 'yyyy-mm-dd HH:MM:SS');
TessMat.ConnEigenmodes = Store;

%% ===== HISTORY + SINGLE SAVE =====
TessMat = bst_history('add', TessMat, 'conn_eigenmodes', ...
    sprintf('Computed %d connection-Laplacian eigenmodes (nSym=%d, reg=%.1e)', ...
        ConnEig.nModes, ConnEig.nSym, ConnEig.Regularization));
bst_save(SurfaceFileFull, TessMat, 'v7');

if isInteractive
    fprintf('BST> Saved %d connection eigenmodes to: %s\n', ConnEig.nModes, SurfaceFile);
end
end
