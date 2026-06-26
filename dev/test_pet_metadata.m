function test_pet_metadata(doEndToEnd)
% TEST_PET_METADATA: validate PET metadata capture on the first 2 PREVENT-AD subjects.
%
% Read-only by default:
%   1. Unit test  - pet_read_metadata on the 2 subjects' BIDS JSON sidecars
%                   (both tracers): assert tracer, units, frame count, coverage.
%   2. Fallback   - pet_read_metadata with no sidecar + a synthetic cube:
%                   assert Source='nifti-header' and frame count from the cube.
%
% Optional (DB-mutating, pass doEndToEnd=1):
%   3. End-to-end - delete existing PET bases for the 2 subjects in the
%                   'preventad' protocol and re-run preventad_pet_import; assert
%                   each surviving base node carries a populated .PET struct.
%
% USAGE:  test_pet_metadata           % read-only checks (safe)
%         test_pet_metadata(1)        % + end-to-end re-import (mutates DB)
%
% Author: Diellor Basha, 2026

    if (nargin < 1) || isempty(doEndToEnd), doEndToEnd = 0; end

    PetRoot  = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet';
    Subjects = {'sub-MTL0002', 'sub-MTL0005'};
    % Expected per tracer: {trc-entity, TracerName, Units, N, [t0 t1] min p.i.}
    Expect = { ...
        '18FNAV4694',      'NAV4694',      'nCi/ml', 6, [40 70]; ...
        '18Fflortaucipir', 'Flortaucipir', 'nCi/ml', 4, [80 100] };

    nPass = 0; nFail = 0;

    %% ===== 1) UNIT TEST: JSON sidecar parsing =====
    fprintf('\n=== PET metadata unit test (BIDS JSON) ===\n');
    for iS = 1:numel(Subjects)
        for iT = 1:size(Expect, 1)
            trc   = Expect{iT,1};
            petFile = fullfile(PetRoot, Subjects{iS}, 'ses-01', 'pet', ...
                               sprintf('%s_ses-01_trc-%s_pet.nii.gz', Subjects{iS}, trc));
            tag = sprintf('%s / %s', Subjects{iS}, trc);
            if ~exist(petFile, 'file')
                [nPass,nFail] = report(nPass, nFail, false, [tag ': PET file present'], 'file missing');
                continue;
            end
            PET = pet_read_metadata(petFile);   % JSON found alongside -> no sMri needed
            [nPass,nFail] = report(nPass,nFail, strcmp(PET.Source,'bids-json'),       [tag ' Source'], PET.Source);
            [nPass,nFail] = report(nPass,nFail, strcmp(PET.Tracer.Name, Expect{iT,2}),[tag ' Tracer'], PET.Tracer.Name);
            [nPass,nFail] = report(nPass,nFail, strcmp(PET.Tracer.Units, Expect{iT,3}),[tag ' Units'],  PET.Tracer.Units);
            [nPass,nFail] = report(nPass,nFail, isequal(PET.Frames.N, Expect{iT,4}),   [tag ' N'],      num2str(PET.Frames.N));
            covOk = ~isempty(PET.Frames.CoverageMinPI) && ...
                    all(abs(PET.Frames.CoverageMinPI - Expect{iT,5}) < 0.1);
            [nPass,nFail] = report(nPass,nFail, covOk, [tag ' Coverage'], mat2str(PET.Frames.CoverageMinPI,4));
            midOk = (numel(PET.Frames.MidTimes) == Expect{iT,4});
            [nPass,nFail] = report(nPass,nFail, midOk, [tag ' MidTimes length'], num2str(numel(PET.Frames.MidTimes)));
            [nPass,nFail] = report(nPass,nFail, ~isempty(PET.Json), [tag ' raw Json kept'], class(PET.Json));
        end
    end

    %% ===== 2) FALLBACK TEST: no sidecar, derive from cube =====
    fprintf('\n=== PET metadata fallback test (no JSON) ===\n');
    fakeFile = fullfile(tempdir, 'no_such_pet.nii.gz');   % guaranteed no .json sidecar
    sMriFake = struct('Cube', zeros(2,2,2,5), 'Header', []);
    PETf = pet_read_metadata(fakeFile, sMriFake);
    [nPass,nFail] = report(nPass,nFail, strcmp(PETf.Source,'nifti-header'), 'fallback Source', PETf.Source);
    [nPass,nFail] = report(nPass,nFail, isequal(PETf.Frames.N, 5),          'fallback N from cube', num2str(PETf.Frames.N));
    [nPass,nFail] = report(nPass,nFail, isempty(PETf.Tracer.Name),          'fallback tracer empty', PETf.Tracer.Name);
    [nPass,nFail] = report(nPass,nFail, isempty(PETf.Json),                 'fallback Json empty', class(PETf.Json));

    %% ===== 3) END-TO-END (optional, mutates DB) =====
    if doEndToEnd
        fprintf('\n=== PET metadata end-to-end re-import ===\n');
        iProtocol = bst_get('Protocol', 'preventad');
        if isempty(iProtocol)
            fprintf('  SKIP: protocol "preventad" not found.\n');
        else
            if (iProtocol ~= bst_get('iProtocol'))
                gui_brainstorm('SetCurrentProtocol', iProtocol);
            end
            for iS = 1:numel(Subjects)
                % delete existing PET base nodes ("PET <tracer>") for a clean re-import
                [sSubject, iSubject] = bst_get('Subject', Subjects{iS});
                if isempty(iSubject), continue; end
                isPetBase = ~cellfun(@isempty, regexp({sSubject.Anatomy.Comment}, '^PET ', 'once'));
                delFiles = {sSubject.Anatomy(isPetBase).FileName};
                if ~isempty(delFiles)
                    file_delete(cellfun(@file_fullpath, delFiles, 'UniformOutput', 0), 1);
                    db_reload_subjects(iSubject);
                end
                % re-import (metadata-aware)
                Out = preventad_pet_import(PetRoot, Subjects{iS}); %#ok<NASGU>
                % verify each surviving base carries .PET
                sSubject = bst_get('Subject', Subjects{iS});
                iBase = find(~cellfun(@isempty, regexp({sSubject.Anatomy.Comment}, '^PET ', 'once')));
                for k = iBase(:)'
                    w = load(file_fullpath(sSubject.Anatomy(k).FileName), 'PET');
                    ok = isfield(w,'PET') && ~isempty(w.PET) && isfield(w.PET,'Source') && strcmp(w.PET.Source,'bids-json');
                    [nPass,nFail] = report(nPass,nFail, ok, ...
                        sprintf('%s / %s carries .PET', Subjects{iS}, sSubject.Anatomy(k).Comment), ...
                        iif(ok, w.PET.Tracer.Name, 'missing'));
                end
            end
        end
    end

    %% ===== summary =====
    fprintf('\n=== SUMMARY: %d passed, %d failed ===\n', nPass, nFail);
    if nFail > 0
        error('test_pet_metadata: %d checks FAILED.', nFail);
    end
end


%% ===== helpers =====
function [nPass, nFail] = report(nPass, nFail, cond, name, got)
    if cond
        nPass = nPass + 1;
        fprintf('  PASS  %s\n', name);
    else
        nFail = nFail + 1;
        fprintf('  FAIL  %s  (got: %s)\n', name, got);
    end
end

function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end
