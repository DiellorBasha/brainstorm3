%use flair instead of T2 - %m51
for k=16:height(sessionJoin)
  
  sourceSubject = sessionJoin.sub{k};
    sourceSession = sessionJoin.MRIsession{k}; %find MRI session
    sourceFolder = inputMRIdir;

    % Main anatomical MRI T1weighted
    inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'T1w');
    % T1w = inputVolume{1};

    % % Secondary MRI FLAIR - for refinement of pial surface
    % inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'FLAIR');
    % 
    if isempty(inputVolume)
        sourceSession(end)='B';
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'T1w');
        if ~isempty (inputVolume)
        T1w = inputVolume{1};
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'FLAIR');
        FLAIR = inputVolume{1};
        end
    else 
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'T1w');
        T1w = inputVolume{1};
        
        FlairInputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'FLAIR');

    end
   
    % % Secondary MRI FLAIR - for refinement of pial surface
        if isempty(FlairInputVolume)
        sourceSession(end)='B';
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'T1w');
       
        T1w = inputVolume{1};
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'FLAIR');
        FLAIR = inputVolume{1};
        end
   
    fs_recon_all (T1w, FLAIR, sourceSubject);
    fs_segmentation(subjectsDir, sourceSubject); %refined segmentation for partial volume corrections

    % get pet data
    sourceSession = sessionJoin.PET_session{k};
    sourceFolder = inputSUVRdir;
    inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'suvr');

    % make templates out of pet SUVRs 
    for jj = 1:length(inputVolume)
        suvr = inputVolume{jj};
        fs_mri_convert (subjectsDir, sourceSubject, suvr);
    end

    a = dir(fullfile(subjectsDir, sourceSubject, 'pet'));
    for jjj = 3:5
        selectedTemplateDir = fullfile(a(jjj).folder, a(jjj).name);

        fs_pet_coreg (subjectsDir, sourceSubject, selectedTemplateDir);
        fs_pet_PVC (subjectsDir, sourceSubject, selectedTemplateDir);
        fs_project_pet (subjectsDir, sourceSubject, selectedTemplateDir);
    end

end

%%
for k=1:25
  
  sourceSubject = sessionJoin.sub{k};
% sourceSubject = 'sub-MTL0005';

    a = dir(fullfile(subjectsDir, sourceSubject, 'pet'));
    for jjj = 3:5
        selectedTemplateDir = fullfile(a(jjj).folder, a(jjj).name);

       % fs_pet_coreg (subjectsDir, sourceSubject, selectedTemplateDir);
        %fs_pet_PVC (subjectsDir, sourceSubject, selectedTemplateDir);
        fs_project_pet (subjectsDir, sourceSubject, selectedTemplateDir);
    end
end
%% % Define the list of subjects
numWorkers = 8; % Adjust based on your system's capabilities
% Set up the parallel pool
parpool(numWorkers);
%% 

% Run gtmseg in parallel

parfor k=2
%fs_T2_pial_refinement(sourceSubject);
sourceSubject = sessionJoin.sub{k};
fs_segmentation(subjectsDir, sourceSubject);
end

% Shut down the parallel pool
delete(gcp('nocreate'));
%% create pet folder, find pet SUVR niftis and make template
k=2
sourceSubject = sessionJoin.sub{k};
sourceSession = sessionJoin.PET_session{k};
sourceFolder = inputSUVRdir;
inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'suvr');

for jj = 1:length(inputVolume)
suvr = inputVolume{jj};

fs_mri_convert (subjectsDir, sourceSubject, suvr);

end

%% co-register

a = dir(fullfile(subjectsDir, sourceSubject, 'pet'))

for k=2
    sourceSubject = sessionJoin.sub{k};

    for jj = 3:5
        selectedTemplateDir = fullfile(a(jj).folder, a(jj).name);
        fs_pet_coreg (subjectsDir, sourceSubject, selectedTemplateDir);
    end
end


%% partial volume correction

a = dir(fullfile(subjectsDir, sourceSubject, 'pet'));

for k=2
    sourceSubject = sessionJoin.sub{k};

    for jj = 3:5
        selectedTemplateDir = fullfile(a(jj).folder, a(jj).name);
        fs_pet_PVC (subjectsDir, sourceSubject, selectedTemplateDir);
    end
end

%% projectPET

a = dir(fullfile(subjectsDir, sourceSubject, 'pet'));

for k=2
    sourceSubject = sessionJoin.sub{k};

    for jj = 3:5
        selectedTemplateDir = fullfile(a(jj).folder, a(jj).name);
        fs_project_pet (subjectsDir, sourceSubject, selectedTemplateDir);
    end
end

%% view template

a = dir(fullfile(subjectsDir, sourceSubject, 'pet'))
selectedTemplate = a(3).name;
freeviewPETtemplate(subjectsDir, sourceSubject, petTemplate);

freeviewRecon(subjectsDir, sourceSubject);

%%
for k=1:height(sessionJoin);
fs_pet_coreg (sourceSubject);
fs_pet_PVC (sourceSubject);
fs_project_pet(sourceSubject);
end
