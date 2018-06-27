function [] = save_all_code(codePathDir, saveFolder)
% Save a copy of all folders and subfolders contained in codePathDir into 
% saveFolder
%
% Notes:
% (1) saveFolder should be the path where everything else about the
%     simulation is being saved

% create the saveFolder if it does not already exist
if ~exist(saveFolder,'dir')
    mkdir(saveFolder); 
end

% make the full path for saving the code folder
savePath = fullfile(saveFolder, 'copy_of_code');

% create the saveFolder if it does not already exist
if ~exist(savePath,'dir')
    mkdir(savePath); 
end

% save a copy of the code
copyfile(codePathDir, savePath);

end

