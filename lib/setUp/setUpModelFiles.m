function [] = setUpModelFiles(modelName,deleteCurrentFiles)
% modelName - a string containitg the model name of the model you want to get
% deleteCurrentFiles - a boolean thet determines if all current model files should be cleared 

if nargin ==2 && deleteCurrentFiles
    % delete all existing model files in Models folder
    delete('./Models/symsFiles/*.m')
    delete('./Models/mex_and_simFiles/*.m')
    delete('./Models/mex_and_simFiles/*.mex*')
end

%find all model symfiles from the development time line
modelsymFiles = dir(sprintf('../developmentTimeline/*%s*/*%s*_syms.m',modelName,modelName));

if ~exist('./Models/symsFiles/','dir') % check if a "symsFiles" directory exists 
   mkdir('./Models/symsFiles/') % if it does not exist create it.
end

% copy over all model files with the correct modelName
for i = 1:numel(modelsymFiles)
    copyfile([modelsymFiles(i).folder,filesep,modelsymFiles(i).name],['./Models/symsFiles/',modelsymFiles(i).name])
end
end