
% Initialize the enviornment for input data and Freesurfer outputs
pipelineName = 'suvr_freesurfer'

fs_initialize
    referenceRegions = '7 8 46 47'; %whole cerebellum - for gray matter do 8 47
    smoothFactor = '5';

for k=36:height(sessionJoin)

sourceSubject = sessionJoin.sub{k};
MRISession = sessionJoin.MRIsession{k};
PETSession = sessionJoin.PET_session{k};

if ~exist(fullfile(subjectsDir, sourceSubject, outputSURFDir), 'dir')
    % do recon if folder does not exist

    T1Files = fs_findInputVolume(sourceSubject, MRISession, '_T1w', 'SourceSubmodality', 'anat', 'BaseDir', inputMRIDir);
    FlairFiles = fs_findInputVolume(sourceSubject, MRISession, '_FLAIR', 'SourceSubmodality', 'anat', 'BaseDir', inputMRIDir);

    T1w = T1Files(1).path;
    FLAIR = FlairFiles(1).path;
    
    if ~isfile(T1w)
        MRISession(end)='B';
            T1Files = fs_findInputVolume(sourceSubject, MRISession, '_T1w', 'SourceSubmodality', 'anat', 'BaseDir', inputMRIDir);
            FlairFiles = fs_findInputVolume(sourceSubject, MRISession, '_FLAIR', 'SourceSubmodality', 'anat', 'BaseDir', inputMRIDir);
                T1w = T1Files(1).path;
                FLAIR = FlairFiles(1).path;
    end
    
        fs_recon_all (T1w, FLAIR, sourceSubject);
        fs_segmentation(subjectsDir, sourceSubject);
end

 % PET Processing
     % PET Processing
 sourceSession = sessionJoin.PET_session{k};
 PETFiles = fs_findInputVolume(sourceSubject, PETSession, '_pet', 'SourceSubmodality', 'pet', 'BaseDir', inputPETDir);
    
   for jj=1:length(PETFiles)
    
    inputPET = PETFiles(jj).path;
   
         inputPETInfo = niftiinfo(inputPET);
         nFrames = inputPETInfo.ImageSize(end);
                
        [submodalityDir, templateDir] = fs_create_templateDir (sourceSubject, inputPET);
      
        disp ("============================================================== \n")
        disp (sprintf ("Processing %s\n", submodalityDir));
        disp ("============================================================== \n")

        registrationFiles = cell(1,nFrames); registeredVolumes = cell(1,nFrames);
       
         for ii = 1:nFrames
            %1. Convert each PET frame to .mgz files
            convertedFrame = fs_frame_convert (templateDir, inputPET, ii-1);
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
                % use brain mask for coregistration
            
            meanPET_regFile = fs_pet_coreg (subjectsDir, sourceSubject, templateDir, meanPET_smoothed);
                % resample to Freesurfer space
            meanPET_resampled = fs_vol2vol (subjectsDir,sourceSubject,templateDir,meanPET_smoothed, meanPET_regFile);

             %====================== PVC ===========
       
       % 5. Mask extracortical regions + rescale to whole cerebellum
       
        PVCFolder = fs_pet_PVC (subjectsDir, sourceSubject, templateDir, meanPET, meanPET_regFile);        
        PVCVolume = 'mgx.ctxgm.nii.gz'; PVCReg = fullfile(PVCFolder, 'aux/bbpet2anat.lta');
        
        % 6. Resample the average to Freesurfer space using mri_vol2vol
        registeredPVC = fs_vol2vol (subjectsDir,sourceSubject,PVCFolder,PVCVolume, PVCReg);

        % Project to surface using mri_vol2surf
        vol2project = fullfile(templateDir, PVCVolume);
        volRegFile = PVCReg;
        [~, ~] = fs_project_pet (subjectsDir, sourceSubject,templateDir, vol2project, volRegFile, PVCFolder);
            
            % ==================== Calculate SUVR without PVC =============

        % 5. Convert averaged frames to SUVR
        
        inputProjVolume = 'template_mean_smoothed_coreg.mgz';
            registrationFile = meanPET_regFile;
            regions = '7 8 46 47'; %whole cerebellum - for gray matter do 8 47

      meanPET_resampled = fs_vol2vol (subjectsDir,sourceSubject,templateDir,meanPET_smoothed, meanPET_regFile);


            suvr = fs_calculate_suvr(subjectsDir, sourceSubject, templateDir, inputProjVolume, registrationFile, regions)
             vol2project = suvr;
             volRegFile = registrationFile;
        
        [r_surf_suvr, l_surf_suvr] = fs_project_pet (subjectsDir, sourceSubject,templateDir, vol2project, volRegFile, templateDir);

        rh_surfs{k} = r_surf_suvr;
        lh_surfs{k} = l_surf_suvr;

   

   end

end