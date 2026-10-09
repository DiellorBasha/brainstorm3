function rec = bst_pet_tracer(iSubject, SubjectName, tracer, method, fwhm, keep4D, OutputDir, bstCommit)
% BST_PET_TRACER  One tracer of the nsp bst-pet worker: static mean, PVC, SUVR, cortex map, regional CSV.
%
% USAGE:  rec = bst_pet_tracer(iSubject, SubjectName, tracer, method, fwhm, keep4D, OutputDir, bstCommit)
%
% Expects the subject to hold the registered base "PET <tracer>" (preventad_pet_import), its T1, and
% the ASEG + Desikan-Killiany volume atlases (import_anatomy_fs). method: 'vlpp' | 'gtm' | 'mg'.
% fwhm: 'vlpp' smoothing in mm; 'gtm'/'mg' scanner FWHM override in mm, or [] for the scanner PSF
% stored with the PET at import (pet_psf_fwhm adds recorded smoothing in quadrature).
% rec: .Tracer .Method .RefValue .GlobalCorticalSUVR .SurfaceFile .RoiCsv .GtmCsv .Provenance
%      (.Method .BrainstormCommit .PsfFwhmMm .PsfSource .SmoothFwhmMm .Reference).
%
% SEE ALSO: bst_pet_subject, pet_pvc, pet_gtm, pet_psf_fwhm, pet_suvr
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-pet pathway)
    rec = struct('Tracer', tracer, 'Method', method, 'RefValue', NaN, ...
                 'GlobalCorticalSUVR', NaN, 'SurfaceFile', '', 'RoiCsv', '', 'GtmCsv', '', 'Provenance', []);
    af = @(c) local_anat_file(iSubject, c);
    baseFile = af(['PET ' tracer]);
    assert(~isempty(baseFile), 'PET %s: registered base missing after import', tracer);
    T1File   = local_t1_file(iSubject);
    sAseg    = in_mri_bst(af('ASEG'));

    % static (mean over frames) — kept as its own node
    sMean = mri_aggregate(in_mri_bst(baseFile), 'mean');
    sMean.Comment = ['PET ' tracer '_mean'];
    meanFile = db_add(iSubject, sMean);

    prov = struct('Method', method, 'BrainstormCommit', bstCommit, 'PsfFwhmMm', [], 'PsfSource', '', ...
                  'SmoothFwhmMm', [], 'Reference', []);
    switch method
        case 'vlpp'   % smooth + plain cerebellar-cortex reference (validated vs VLPP)
            sIn = sMean;
            sIn.Cube = imgaussfilt3(double(sMean.Cube(:,:,:,1)), fwhm / 2.355);
            suvrOpts = struct('Erode', 0, 'Robust', 'mean');
            prov.SmoothFwhmMm = fwhm;
        case {'gtm', 'mg'}   % PVC (pr/pet-pvc-2 API), robust eroded reference
            PET = [];
            if isfield(sMean, 'PET'), PET = sMean.PET; end
            % The PSF pet_pvc resolves itself from the same PET metadata: recorded here for provenance
            [prov.PsfFwhmMm, prov.PsfSource] = pet_psf_fwhm(PET, fwhm);
            if isstruct(PET) && isfield(PET, 'SmoothFwhm'), prov.SmoothFwhmMm = PET.SmoothFwhm; end
            [pvcFile, errPvc, ~, regTable] = pet_pvc(meanFile, T1File, fwhm, struct('method', method));
            assert(isempty(errPvc) && ~isempty(pvcFile), 'PET %s: pet_pvc %s failed: %s', tracer, method, errPvc);
            sIn = in_mri_bst(pvcFile);  suvrOpts = struct();
            if strcmp(method, 'gtm')
                rec.GtmCsv = local_write_gtm(regTable, fullfile(OutputDir, ...
                    sprintf('%s_pet_gtm_%s.csv', SubjectName, tracer)));
            end
    end
    [sSuvr, info] = pet_suvr(sIn, sAseg, suvrOpts);
    prov.Reference = struct('Region', 'cerebellar cortex', 'AsegLabels', [8 47], 'Mode', info.Reference, ...
                            'Erode', info.Erode, 'Estimator', info.Robust, 'Value', info.RefValue, ...
                            'nVoxMask', info.nVoxMask, 'nVoxUsed', info.nVoxEroded);
    rec.Provenance = prov;
    sSuvr.Comment = ['PET ' tracer '_suvr'];
    sSuvr = bst_history('add', sSuvr, 'suvr', sprintf( ...
        'SUVR (%s), ref cerebellar cortex = %.4g', method, info.RefValue));
    suvrFile = db_add(iSubject, sSuvr);
    rec.RefValue = info.RefValue;

    % cortical map on the subject's cortex (condition "PET")
    [surfFile, errProj] = mri_interp_vol2tess(suvrFile, T1File, 'PET', 'SUVR', [0.1 0.8 0.1]);
    if ~isempty(errProj), warning('PET %s: projection failed: %s', tracer, errProj); end
    rec.SurfaceFile = surfFile;

    % regional (Desikan) + global cortical SUVR
    sDK = in_mri_bst(af('Desikan-Killiany'));
    [names, vals] = local_regional(sDK, double(sSuvr.Cube(:,:,:,1)));
    rec.GlobalCorticalSUVR = mean(vals, 'omitnan');
    rec.RoiCsv = fullfile(OutputDir, sprintf('%s_pet_suvr_%s.csv', SubjectName, tracer));
    fid = fopen(rec.RoiCsv, 'w');
    fprintf(fid, 'subject,tracer,method,reference,global_cortical_suvr');
    fprintf(fid, ',%s', names{:}); fprintf(fid, '\n');
    fprintf(fid, '%s,%s,%s,%.6g,%.6g', SubjectName, tracer, method, info.RefValue, rec.GlobalCorticalSUVR);
    fprintf(fid, ',%.6g', vals); fprintf(fid, '\n');
    fclose(fid);
    [~, n, e] = fileparts(rec.RoiCsv); rec.RoiCsv = [n e];

    if ~keep4D
        local_delete_anat(iSubject, baseFile);
    end
    fprintf('PET %s: ref=%.4g global cortical SUVR=%.3f -> %s\n', ...
        tracer, info.RefValue, rec.GlobalCorticalSUVR, surfFile);
end


%% ===== helpers =====

function f = local_anat_file(iSubject, comment)
    sSubject = bst_get('Subject', iSubject);
    i = find(strcmp({sSubject.Anatomy.Comment}, comment), 1);
    if isempty(i), f = ''; else, f = sSubject.Anatomy(i).FileName; end
end

function f = local_t1_file(iSubject)
    sSubject = bst_get('Subject', iSubject);
    f = sSubject.Anatomy(sSubject.iAnatomy).FileName;
end

function local_delete_anat(iSubject, fileName)
    sSubject = bst_get('Subject', iSubject);
    i = find(strcmp({sSubject.Anatomy.FileName}, fileName), 1);
    if isempty(i), return; end
    file_delete(file_fullpath(fileName), 1);
    sSubject.Anatomy(i) = [];
    if sSubject.iAnatomy > i, sSubject.iAnatomy = sSubject.iAnatomy - 1; end
    bst_set('Subject', iSubject, sSubject);
end

function [names, vals] = local_regional(sDK, cube)
% Desikan-Killiany cortical labels (1000-2999), as in preventad_pet_pipeline.
    L = sDK.Labels; v = cell2mat(L(:,1)); nm = L(:,2);
    k = find(v >= 1000 & v < 3000);
    names = nm(k)'; vals = zeros(1, numel(k));
    for i = 1:numel(k), vals(i) = mean(cube(sDK.Cube == v(k(i))), 'omitnan'); end
end

function name = local_write_gtm(T, file)
% pet_gtm regTable -> CSV (id,name,nvox,observed,corrected); returns the file name only.
    fid = fopen(file, 'w');
    fprintf(fid, 'id,name,nvox,observed,corrected\n');
    for i = 1:numel(T)
        fprintf(fid, '%d,"%s",%d,%.6g,%.6g\n', T(i).id, T(i).name, T(i).nvox, T(i).observed, T(i).corrected);
    end
    fclose(fid);
    [~, n, e] = fileparts(file); name = [n e];
end
