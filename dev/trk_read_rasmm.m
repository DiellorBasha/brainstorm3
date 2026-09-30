function [tracks, header] = trk_read_rasmm(FileName)
% TRK_READ_RASMM  Read a TrackVis .trk as TrackVis voxmm, with the header's vox_to_ras.
%
% USAGE:  [tracks, header] = trk_read_rasmm(FileName)
%
% Unlike external/trk/trk_read.m, the points are returned as STORED (TrackVis
% "voxmm": millimetres from the corner of the .trk's voxel grid) with no
% orientation flips, and the header carries the 4x4 vox_to_ras that places them in
% scanner/world RAS millimetres. Convert with trk_voxmm_to_rasmm, after any
% resampling (trk_interp's cubic B-splines are affine-invariant, so resampling in
% voxmm and then transforming is exact).
%
% OUTPUTS:
%    - tracks : struct array [1 x nTracks] with .nPoints and .matrix [nPoints x 3]
%               (the layout trk_interp expects); scalars and properties are dropped
%    - header : .dim, .voxel_size, .vox_to_ras (4x4), .voxel_order (e.g. 'LPS'),
%               .n_scalars, .n_properties, .n_count
%
% The TrackVis v2 header layout (1000 bytes, little-endian): dim int16[3] @6,
% voxel_size float32[3] @12, n_scalars int16 @36, n_properties int16 @238,
% vox_to_ras float32[4x4] row-major @440, voxel_order char[4] @948,
% n_count int32 @988, version int32 @992, hdr_size int32 @996.
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-fibers pathway)

fid = fopen(FileName, 'r', 'ieee-le');
if fid < 0
    error('trk_read_rasmm: cannot open %s', FileName);
end
cleaner = onCleanup(@() fclose(fid));

raw = fread(fid, 1000, '*uint8')';
if numel(raw) < 1000 || ~strcmp(char(raw(1:5)), 'TRACK')
    error('trk_read_rasmm: %s is not a TrackVis file', FileName);
end
at = @(off, n, type) typecast(raw(off+1 : off+n*bytes(type)), type);
header.dim          = double(at(6, 3, 'int16'));
header.voxel_size   = double(at(12, 3, 'single'));
header.n_scalars    = double(at(36, 1, 'int16'));
header.n_properties = double(at(238, 1, 'int16'));
header.vox_to_ras   = reshape(double(at(440, 16, 'single')), 4, 4)';   % stored row-major
header.voxel_order  = upper(strtrim(char(raw(949:951))));
header.n_count      = double(at(988, 1, 'int32'));
header.version      = double(at(992, 1, 'int32'));
if double(at(996, 1, 'int32')) ~= 1000
    error('trk_read_rasmm: %s has a header size other than 1000 bytes', FileName);
end
if header.version ~= 2 || all(header.vox_to_ras(:) == 0)
    error('trk_read_rasmm: %s has no vox_to_ras (TrackVis v2 required)', FileName);
end
% the voxel order must be the affine's own orientation: then no axis reordering applies
affOrder = local_orientation(header.vox_to_ras);
if ~isempty(header.voxel_order) && ~strcmp(header.voxel_order, affOrder)
    error('trk_read_rasmm: voxel_order %s differs from vox_to_ras orientation %s', header.voxel_order, affOrder);
end

nCols = 3 + header.n_scalars;
nMax  = header.n_count;
if nMax <= 0; nMax = Inf; end
tracks = repmat(struct('nPoints', 0, 'matrix', []), 1, min(nMax, 1e7));
iTrk = 0;
while iTrk < nMax
    n = fread(fid, 1, 'int32');
    if isempty(n) || feof(fid); break; end
    iTrk = iTrk + 1;
    m = fread(fid, [nCols, n], '*single')';
    tracks(iTrk).nPoints = n;
    tracks(iTrk).matrix  = m(:, 1:3);
    if header.n_properties
        fread(fid, header.n_properties, 'single');
    end
end
tracks = tracks(1:iTrk);
header.n_count = iTrk;
end

function n = bytes(type)
switch type
    case 'int16';  n = 2;
    case 'int32';  n = 4;
    case 'single'; n = 4;
end
end

function s = local_orientation(A)
% The axis codes of a voxel-to-RAS affine (nibabel's aff2axcodes): per voxel axis,
% the world axis it runs along most and its sign.
labels = ['LR'; 'PA'; 'IS'];
s = '';
R = A(1:3, 1:3);
for j = 1:3
    [~, i] = max(abs(R(:, j)));
    s(j) = labels(i, 1 + (R(i, j) > 0)); %#ok<AGROW>
end
end
