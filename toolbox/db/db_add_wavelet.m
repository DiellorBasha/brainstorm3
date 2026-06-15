function iWavelet = db_add_wavelet(iSubject, ParentEigenFile, WaveletMat, Comment)
% DB_ADD_WAVELET: Save a wavelet_*.mat and register it as a child of an eigen node.
%
% USAGE:  iWavelet = db_add_wavelet(iSubject, ParentEigenFile, WaveletMat, Comment)
%
% INPUT:
%    - iSubject        : Index of the subject (0 = default subject); reconciled against
%                        the subject resolved from the parent eigen node.
%    - ParentEigenFile : Relative or full path to the parent eigen_*.mat node.
%    - WaveletMat      : Structure to save (db_template('waveletmat') schema).
%    - Comment         : Description string for this wavelet node (default: 'Wavelet').
% OUTPUT:
%    - iWavelet        : Index of the new entry in sSubject.Surface(iSurface).Wavelet
%
% The node is stored in the SURFACE's Wavelet list, keyed by ParentEigen (the eigen
% node file). node_create_subject nests it under the matching eigen node in the tree.
%
% Authors: Diellor Basha, 2026

    if (nargin < 4) || isempty(Comment); Comment = 'Wavelet'; end

    % Resolve the parent eigen node -> its subject/surface
    [~, iSubjectE, iSurface] = bst_get('EigenFile', ParentEigenFile);
    if isempty(iSurface)
        error('db_add_wavelet:noEigen', 'Parent eigen node not found: %s', ParentEigenFile);
    end
    iSubject = iSubjectE;   % trust the resolved subject (the argument is kept for API symmetry)

    ProtocolSubjects = bst_get('ProtocolSubjects');
    if (iSubject == 0); sSubject = ProtocolSubjects.DefaultSubject;
    else;               sSubject = ProtocolSubjects.Subject(iSubject); end

    % Build a unique filename in the parent surface's anatomy folder
    ProtocolInfo = bst_get('ProtocolInfo');
    c = clock;
    strTime = sprintf('%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
    OutputFile = ['wavelet_' strTime '.mat'];
    OutputFileFull = file_unique(bst_fullfile(ProtocolInfo.SUBJECTS, bst_fileparts(sSubject.FileName), OutputFile));

    % Stamp required fields
    WaveletMat.ParentEigen = file_short(ParentEigenFile);
    WaveletMat.Comment     = Comment;

    bst_save(OutputFileFull, WaveletMat, 'v7');

    % Normalize the surface array to the current template (adds Wavelet if missing),
    % mirroring db_add_eigen's homogeneity fix (a fresh array is rebuilt because widening
    % the field set of one element of an existing struct array in place throws).
    templateSurface = db_template('Surface');
    tFields = fieldnames(templateSurface);
    if ~isempty(sSubject.Surface) && ~all(isfield(sSubject.Surface, tFields))
        sNew = cell(1, numel(sSubject.Surface));
        for k = 1:numel(sSubject.Surface)
            s = struct_copy_fields(sSubject.Surface(k), templateSurface, 0);
            extraFields = setdiff(fieldnames(s), tFields, 'stable');
            sNew{k} = orderfields(s, [tFields; extraFields]);
        end
        sSubject.Surface = reshape([sNew{:}], size(sSubject.Surface));
    end

    % Append the child entry
    newEntry.FileName    = file_short(OutputFileFull);
    newEntry.Comment     = Comment;
    newEntry.ParentEigen = file_short(ParentEigenFile);
    sSubject.Surface(iSurface).Wavelet(end+1) = newEntry;
    iWavelet = numel(sSubject.Surface(iSurface).Wavelet);

    % Write the modified subject back via the canonical round-trip
    if (iSubject == 0); ProtocolSubjects.DefaultSubject = sSubject;
    else;               ProtocolSubjects.Subject(iSubject) = sSubject; end
    bst_set('ProtocolSubjects', ProtocolSubjects);

    % Refresh tree and save database
    panel_protocols('UpdateNode', 'Subject', iSubject);
    db_save();
end
