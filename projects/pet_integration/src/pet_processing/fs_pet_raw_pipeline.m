inputMRIdir = '/export02/data/dbasha/PAD7/data_release_7.0/mri/wave1'
inputPETdir =  '/export02/data/dbasha/PAD7/data_release_7.0/pet';
inputSUVRdir = '/export02/data/dbasha/PAD7/data_release_7.0/pet/derivatives/vlpp';
stormPath = '/export01/data/dbasha/code/StormNet';
brainstormPath = '/meg/meg1/software/brainstorm3'
    addpath(genpath(stormPath));
    addpath(genpath(brainstormPath));
% ==============
%% 

% Project SUVRs to surface without registration (they are already
% registered)
%
% Process PET to get SUVR then project to surface
%  incr = [1:24]
% [27:height(sessionJoin)

% 1. Convert PET frames to .mgz, align them to T1, average the result

for k=26:height(sessionJoin)
    
    sourceSubject = sessionJoin.sub{k};
    checkFolder = fullfile(subjectsDir, sourceSubject);

    if ~isfolder (checkFolder)
        fs_recon_segment (subjectsDir, sourceSubject, sessionJoin, inputMRIdir, k)
    end

    % PET Processing
    sourceSession = sessionJoin.PET_session{k};
    sourceFolder = inputPETdir;
    inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'pet');

    for jj = 1:length(inputVolume)
 
        currentPET = inputVolume{jj};
        currentPETInfo = niftiinfo(currentPET);
        nFrames = currentPETInfo.ImageSize(end);

        templateDir = fs_create_templateDir (subjectsDir, sourceSubject, currentPET);
        registrationFiles = cell(1,nFrames);
        registeredVolumes = cell(1,nFrames);

        for ii = 1:nFrames

            %1. Convert each PET frame to .mgz files
            convertedFrame = fs_frame_convert (templateDir, currentPET, ii-1);

            % 2. Align each PET frame to T1 MRI using mri_coreg

            registrationFile = fs_pet_coreg (subjectsDir, sourceSubject, templateDir, convertedFrame);
            % resample to T1
            registeredVolume = fs_vol2vol (subjectsDir,sourceSubject,templateDir,convertedFrame, registrationFile);
            % store results
            registrationFiles{ii} = registrationFile;
            registeredVolumes{ii} = registeredVolume;
        end

        % 3. Average the PET frames using mri_concat
        subFrames = registeredVolumes(1:nFrames);
        meanPET = fs_mri_concat_frames (subjectsDir, sourceSubject, templateDir, subFrames);

        % 4. coregister the mean to T1
        
        smoothFactor = '5';
        meanPET_smoothed = fs_smooth(subjectsDir,sourceSubject, templateDir, meanPET, smoothFactor);
        meanPET_smoothed2 = 'template_mean_smoothed.mgz';
        % use brain mask for coregistration
        meanPET_regFile = fs_pet_coreg (subjectsDir, sourceSubject, templateDir, meanPET_smoothed2);
        meanPET_resampled = fs_vol2vol (subjectsDir,sourceSubject,templateDir,meanPET_smoothed2, meanPET_regFile);

        %====================== PVC ===========
        % % 5. Masking extracortical regions + rescaling to whole cerebellum - no PVC
        PVCFolder = fs_pet_PVC (subjectsDir, sourceSubject, templateDir, meanPET, meanPET_regFile);
        %
        % % 6. Resample the average to Freesurfer space using mri_vol2vol
        PVCVolume = 'mgx.ctxgm.nii.gz'; PVCReg = fullfile(PVCFolder, 'aux/bbpet2anat.lta');
        %
        registeredVolume = fs_vol2vol (subjectsDir,sourceSubject,PVCFolder,PVCVolume, PVCReg)

               % Project to surface using mri_vol2surf
        vol2project = fullfile(templateDir, PVCVolume);
        volRegFile = PVCReg;
        [~, ~] = fs_project_pet (subjectsDir, sourceSubject,templateDir, vol2project, volRegFile, PVCFolder);


        %        PVC is causing edge artefacts;
        % ==================== PVC =============

        % 7. Convert averaged frames to SUVR
        inputProjVolume = 'template_mean_smoothed_coreg.mgz';
        registrationFile = meanPET_regFile;
        regions = '7 8 46 47'; %whole cerebellum - for gray matter do 8 47

        suvr = fs_calculate_suvr(subjectsDir, sourceSubject, templateDir, inputProjVolume, registrationFile, regions)
        vol2project = suvr;
        volRegFile = registrationFile;
        [r_surf_suvr, l_surf_suvr] = fs_project_pet (subjectsDir, sourceSubject,templateDir, vol2project, volRegFile, templateDir);

        rh_surfs{k} = r_surf_suvr;
        lh_surfs{k} = l_surf_suvr;
    end
        


        % %======== Project VLPP SUVRs for the same subject =========================
    sourceFolder = inputSUVRdir;
    inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'suvr');

        for jj = 1:length(inputVolume)
            currentPET = inputVolume{jj};
            currentPETInfo = niftiinfo(currentPET);

            outputSUVRDir = extractBetween(currentPET, '_trc-', '.nii.gz');
                     
           templateDir = fs_create_templateDir (subjectsDir, sourceSubject, outputSUVRDir{1});


            %1. Convert each PET frame to .mgz files

            convertedVolume = fullfile(templateDir, 'suvr.mgz');
            cmd = ['mri_convert ' currentPET ' ' convertedVolume];
            fs_execute(cmd)

            % 2. Align each PET frame to T1 MRI using mri_coreg

            registrationFile = fs_pet_coreg (subjectsDir, sourceSubject, templateDir, 'suvr.mgz');
            registeredVolume = fs_vol2vol (subjectsDir,sourceSubject,templateDir,'suvr.mgz', registrationFile);

            % Project to surface using mri_vol2surf
            vol2project = registeredVolume;
            volRegFile = registrationFile;
            [r_surf_suvr, l_surf_suvr] = fs_project_pet (subjectsDir, sourceSubject,templateDir, vol2project, volRegFile, templateDir);

            rh_surfs_vlpp{k} = r_surf_suvr;
            lh_surfs_vlpp{k} = l_surf_suvr;
cd(subjectsDir)
        end
end

    
%%

k=2
tt = 2
selectedPETVol = 'suvr.mgz';
selectedPVCVol = 'suvr_gtmpvc/mgx.ctxgm_coreg.mgz';

    sourceSubject = sessionJoin.sub{k};
    tracers = {'18Fflortaucipir_pet'; '18FNAV4694_pet'; '18FNAV4694_pet_ref-cerebellumCortex_suvr'};
  %  templateDir = fullfile(subjectsDir, 'pet', tracers{tt});
    partialTemplateDir = fullfile ('pet', tracers{tt});

avPVC = strcat('./', partialTemplateDir, '/', 'mgx.ctxgm_coreg.mgz:colormap=heat:smoothed=1:heatscale=0,2.5,5')
avPET = strcat('./', partialTemplateDir, '/', 'template_mean_suvr.mgz')
l_overlay = strcat('./', partialTemplateDir, '/lh_surf_suvr.mgz');
r_overlay = strcat('./', partialTemplateDir, '/rh_surf_suvr.mgz');
overlaySett = [':overlay_method=linear' ...
    ':overlay_color=heat,truncate' ...
    ':overlay_threshold=1,3' ...
    ':edgecolor=overlay' ...
    ':edgethickness=0'];
PVCSett = [':colormap=heat' ...
    ':smoothed=1' ...
    ':heatscale=0,5,10' ...
    ':heatscale_options=truncate' ...
    ':mask=brainmask'];
screenshot = ' -screenshot %name 2 auto_trim';

cd (fullfile(subjectsDir, sourceSubject))
cmd = ['freeview' ...
    ' -volume ./mri/brainmask.mgz' ...
    ' -volume ./mri/T1.mgz'...
    ' -volume ' avPVC PVCSett ...
    ' -volume ' avPET PVCSett ...
    ' -surface ./surf/lh.pial:overlay=' l_overlay overlaySett...
    ' -surface ./surf/rh.pial:overlay=' r_overlay overlaySett...
    ' -view right' ...
    ' -hide-3d-slices' ...
    ' -zoom 1.5' ...
    ' -colorscale'];

fs_execute(cmd)

cd(subjectsDir)
