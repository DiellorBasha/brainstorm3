%use flair instead of T2 - %m51 78
function fs_recon_segment (sourceSubject, inputT1, inputFLAIR)
    
    if isempty(inputT1)
        sourceSession(end)='B';
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'T1w');
        if ~isempty (inputVolume)
        T1w = inputVolume{1};
        FlairInputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'FLAIR');
        FLAIR = FlairInputVolume{1};
        end
    else 
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'T1w');
        T1w = inputVolume{1};
        
        FlairInputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'FLAIR');
        FLAIR = FlairInputVolume{1};
    end
   
    % % Secondary MRI FLAIR - for refinement of pial surface
        if isempty(FlairInputVolume)
        sourceSession(end)='B';
        inputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'T1w');
       
        T1w = inputVolume{1};
        FlairInputVolume = findInputVolume (sourceFolder, sourceSubject, sourceSession, 'FLAIR');
        FLAIR = FlairInputVolume{1};
        end
   
    fs_recon_all (T1w, FLAIR, sourceSubject);
    fs_segmentation(subjectsDir, sourceSubject);
end