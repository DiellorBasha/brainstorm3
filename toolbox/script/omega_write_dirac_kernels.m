function OutFiles = omega_write_dirac_kernels()
% OMEGA_WRITE_DIRAC_KERNELS  Write the module's Dirac-eigenbasis amplitude currentKernel back into
% the omega_tutorial_test protocol as a standard unconstrained ImagingKernel result, per subject.
%
% The nxr-cortical-flow-matlab module builds the relative-Dirac eigenbasis inverse and reconstructs
% the vertex current operator currentKernel = reconstruct(diracInverse, dbasis) [3V x C]. That is a
% drop-in for the vanilla wMNE ImagingKernel (same 270 MEG channels, same units, nComponents=3), so
% process_source_flow / process_source_flow_functional fuse the flow operators onto it unchanged --
% giving the EIGENBASIS-regularized flow instead of the plain-wMNE flow.
%
% Clones the existing MN kernel as a template (guarantees all required fields + channel order) and
% swaps ImagingKernel. Channel order was verified identical (ctx.iSel == GoodChannel, sorted).
%
% Author: Diellor Basha, 2026

    addpath('/Users/diellorbasha/workspace/research/code/nxr-cortical-flow-matlab');
    gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', 'omega_tutorial_test'));

    % subject -> (module dataset name, its raw study, the MN-kernel template file)
    spec = {
        'sub-0002', 'omega_sub0002', '@rawsub-0002_ses-01_task-rest_run-01_meg_notch_high', 'results_MN_MEG_KERNEL_260728_2139.mat'
        'sub-0003', 'omega_sub0003', '@rawsub-0003_ses-01_task-rest_run-01_meg_notch_high', 'results_MN_MEG_KERNEL_260728_2139.mat'
    };
    OutFiles = cell(size(spec,1),1);
    touchedStudies = [];

    for i = 1:size(spec,1)
        subj = spec{i,1};  ds = spec{i,2};  raw = spec{i,3};  mnName = spec{i,4};
        mnFile = bst_fullfile(subj, raw, mnName);

        ctx = flow.context(ds, 400);                       % amplitude Dirac currentKernel [3V x C]
        R   = in_bst_results(mnFile, 0);                    % template (all fields, right channel order)
        assert(isequal(ctx.iSel(:), double(R.GoodChannel(:))), 'channel order mismatch for %s', subj);
        assert(isequal(size(ctx.currentKernel), size(R.ImagingKernel)), 'kernel shape mismatch for %s', subj);

        R.ImagingKernel = ctx.currentKernel;                % the swap
        R.Comment       = 'Dirac-flow: MEG (Unconstr) amplitude';
        R.Function      = 'mn';
        R = bst_history('add', R, 'compute', ...
            'Dirac-eigenbasis amplitude inverse (nxr-cortical-flow-matlab source.dirac), currentKernel = reconstruct(diracInverse, dbasis)');

        [~, iStudy] = bst_get('AnyFile', mnFile);
        kdir    = bst_fileparts(file_fullpath(mnFile));
        NewFile = file_unique(bst_fullfile(kdir, 'results_DiracFlowAmp_KERNEL.mat'));
        bst_save(NewFile, R, 'v6');
        db_add_data(iStudy, file_short(NewFile), R);
        OutFiles{i} = file_short(NewFile);
        touchedStudies(end+1) = iStudy; %#ok<AGROW>
        fprintf('[%s] wrote Dirac-flow kernel -> %s\n', subj, OutFiles{i});
    end

    db_reload_studies(unique(touchedStudies));
    panel_protocols('UpdateTree');
    fprintf('omega_write_dirac_kernels: %d kernels written.\n', numel(OutFiles));
end
