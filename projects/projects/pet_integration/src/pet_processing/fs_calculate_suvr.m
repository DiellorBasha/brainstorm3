function suvr = fs_calculate_suvr(subjectsDir, sourceSubject, templateDir, inputVolume, registrationFile, regions)
     
% Get cerebellar suv from GTM stats and do ratio from template_mean
      
% get gtmseg into 256x256x256
        mriDir = fullfile(subjectsDir, sourceSubject, 'mri');
        gtmFile = 'gtmseg.mgz';
        gtm_resampled = fs_vol2vol(subjectsDir, sourceSubject, mriDir, gtmFile,registrationFile);
  
      
 segStatInputVol = fullfile(templateDir, inputVolume);
 segStatSeg = fullfile(subjectsDir, sourceSubject, '/mri/gtmseg_coreg.mgz');
 segStatOutput = fullfile(templateDir, 'mean_uptake_cerebellum.txt')
      
 % Get average in region      
 avg = fs_segstat_regions(segStatInputVol, segStatSeg, segStatOutput, regions);
       
 % Calculate SUVR
 outputSUVRFile = fullfile(templateDir, 'template_mean_suvr.mgz');
 avgString = num2str(avg);
        
        cmd = ['mris_calc' ... 
               ' --output ' outputSUVRFile ...
               ' ' segStatInputVol ...
               ' div ' avgString ...
               ];
        
        fs_execute(cmd)

        suvr = outputSUVRFile; 
end