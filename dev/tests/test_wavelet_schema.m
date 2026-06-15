function test_wavelet_schema()
% Schema regression for the waveletmat template and the Surface.Wavelet list.
% Authors: Diellor Basha, 2026
    nFail = 0;
    w = db_template('waveletmat');
    need = {'Comment','ParentEigen','Variant','Tiles','Tiling','Provenance'};
    for f = need
        if ~isfield(w, f{1}); fprintf('MISSING waveletmat.%s\n', f{1}); nFail = nFail+1; end
    end
    s = db_template('surface');
    if ~isfield(s, 'Wavelet'); fprintf('MISSING surface.Wavelet\n'); nFail = nFail+1; end
    if isfield(s,'Wavelet')
        sub = fieldnames(s.Wavelet);
        for f = {'FileName','Comment','ParentEigen'}
            if ~ismember(f{1}, sub); fprintf('MISSING surface.Wavelet.%s\n', f{1}); nFail = nFail+1; end
        end
        if ~isempty(s.Wavelet); fprintf('surface.Wavelet must start EMPTY\n'); nFail = nFail+1; end
    end
    if isfield(s, 'Filterbank'); fprintf('surface.Filterbank should be GONE\n'); nFail = nFail+1; end
    fprintf('\n==== test_wavelet_schema: %d failed ====\n', nFail);
    if nFail > 0, error('test_wavelet_schema FAILED'); end
end
