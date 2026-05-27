function [Coeffs, Reconstructed] = tess_eigenmodes_project(SurfaceFile, Data, varargin)
% TESS_EIGENMODES_PROJECT: Project scalar data onto eigenmode basis and optionally reconstruct.
%
% USAGE:  Coeffs = tess_eigenmodes_project(SurfaceFile, Data)
%         Coeffs = tess_eigenmodes_project(SurfaceFile, Data, 'ModeRange', [1 100])
%         [Coeffs, Reconstructed] = tess_eigenmodes_project(SurfaceFile, Data)
%
% DESCRIPTION:
%     Projects one or more scalar fields (defined on mesh vertices) onto the
%     precomputed Laplace-Beltrami eigenmode basis. Each field u is decomposed
%     into spectral coefficients:
%
%         c_k = phi_k' * M * u     (M-weighted inner product)
%
%     where phi_k are the eigenmodes and M is the mass matrix. The coefficients
%     c_k represent the "spatial frequency content" of the field — low indices
%     correspond to smooth, large-scale patterns; high indices to fine detail.
%
%     Optionally reconstructs the field from a subset of modes:
%
%         u_approx = Phi(:, ModeRange) * c(ModeRange)
%
%     This is spatial low-pass, band-pass, or high-pass filtering depending
%     on which ModeRange is selected.
%
% INPUTS:
%     SurfaceFile : Relative or absolute path to the Brainstorm surface file
%                   containing precomputed Eigenmodes (from out_tess_eigenmodes).
%     Data        : [nVertices x nTime] or [nVertices x 1] scalar field(s) on
%                   the mesh. Each column is one time point or sample.
%
% OPTIONS:
%     ModeRange   : [1 x 2] range of mode indices to use, e.g., [1 100].
%                   Default: all available modes.
%     MassMatrix  : [nV x nV] sparse mass matrix. If empty, recomputed from
%                   the surface via tess_laplacian. Default: [].
%
% OUTPUTS:
%     Coeffs       : [nModes x nTime] spectral coefficients c_k(t).
%     Reconstructed: [nVertices x nTime] field reconstructed from the selected
%                    mode range. Only computed if requested (2 output args).
%
% MATHEMATICAL NOTES:
%     - The eigenmodes satisfy M-orthonormality: Phi' * M * Phi = I
%     - Parseval's identity: ||u||^2_M = sum(c_k^2)  (exact when all modes used)
%     - Energy in band [k1, k2]: E_band = sum(c_k^2) for k in [k1, k2]
%     - The eigenvalues lambda_k are spatial frequencies squared. The "spatial
%       frequency" of mode k is sqrt(lambda_k).
%
% EXAMPLES:
%     % Project source activity onto eigenmodes
%     Coeffs = tess_eigenmodes_project(SurfaceFile, SourceData);
%
%     % Reconstruct using only the first 50 modes (spatial low-pass)
%     [~, Smoothed] = tess_eigenmodes_project(SurfaceFile, SourceData, 'ModeRange', [1 50]);
%
%     % Band-pass: modes 20-80 (medium spatial frequencies only)
%     [~, BandPass] = tess_eigenmodes_project(SurfaceFile, SourceData, 'ModeRange', [20 80]);
%
% SEE ALSO: tess_eigenmodes, in_tess_eigenmodes, tess_eigenmodes_filter, tess_laplacian

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

%% ===== PARSE INPUTS =====
ModeRange  = [];
MassMatrix = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'moderange',   ModeRange  = varargin{i+1};
        case 'massmatrix',  MassMatrix = varargin{i+1};
    end
end

%% ===== LOAD EIGENMODES =====
[Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile);
if ~isComputed
    error('No precomputed eigenmodes found in: %s. Run tess_eigenmodes first.', SurfaceFile);
end

Phi     = Eigenmodes.Vectors;   % [nV x nModes], double (converted from single by _load)
lambdas = Eigenmodes.Values;    % [nModes x 1]
nModes  = Eigenmodes.nModes;
nV      = size(Phi, 1);

%% ===== VALIDATE DATA =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but surface has %d vertices.', size(Data, 1), nV);
end

%% ===== MODE RANGE =====
if isempty(ModeRange)
    ModeRange = [1, nModes];
end
ModeRange(1) = max(1, ModeRange(1));
ModeRange(2) = min(nModes, ModeRange(2));
iModes = ModeRange(1):ModeRange(2);

if isempty(iModes)
    error('Empty mode range [%d, %d] (surface has %d modes).', ModeRange(1), ModeRange(2), nModes);
end

%% ===== GET OR COMPUTE MASS MATRIX =====
if isempty(MassMatrix)
    % Load surface for mass matrix computation
    SurfaceFileFull = file_fullpath(SurfaceFile);
    if isempty(SurfaceFileFull)
        SurfaceFileFull = SurfaceFile;
    end
    TessMat = load(SurfaceFileFull, 'Vertices', 'Faces');
    [~, MassMatrix] = tess_laplacian(double(TessMat.Vertices), double(TessMat.Faces), ...
        'MassType', Eigenmodes.MassType);
end

%% ===== PROJECT: c_k = phi_k' * M * u =====
% For all modes first (needed even if ModeRange is a subset)
% M * Data is [nV x nTime], Phi' * (M * Data) is [nModes x nTime]
MData = MassMatrix * Data;
Coeffs = Phi' * MData;  % [nModes x nTime]

%% ===== RECONSTRUCT (if requested) =====
if nargout >= 2
    Reconstructed = Phi(:, iModes) * Coeffs(iModes, :);  % [nV x nTime]
end

% Trim Coeffs to requested range if user only wants those
% (Keep all coefficients by default — more useful for analysis)

end
