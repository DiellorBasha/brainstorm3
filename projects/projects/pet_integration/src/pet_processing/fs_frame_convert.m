function convertedVolume = fs_frame_convert (templateDir, inputPETVolume, frameNumber);       
ii=frameNumber;
%label frames starting from iterator 1:nFrames;

outputVolume = fullfile(templateDir, sprintf('frame_%d.mgz', ii+1));

% Pet surfer iterates from 0:nFrames-1 however, we label 1:nFrames
currentFrame = sprintf(' -nth %d', ii);

                cmd = ['mri_convert '...
                        inputPETVolume ...
                        ' ' ...
                        outputVolume ...
                        currentFrame
                ];

        fs_execute(cmd)

        convertedVolume = sprintf('frame_%d.mgz', ii+1);

end
