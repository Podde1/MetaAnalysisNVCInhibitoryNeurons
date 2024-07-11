function [theta,bestParamFileName] = getBestSavedParameters(resultFolderName,modelName,varargin)
%% load parameters for inital values
% Algorithm will prioritise parameters from:
% 1. ./ResultLocation - lowest fval
% 2. ./Parameters/resultFolderName - lowset fval
% 3. ./Parameters - lowest fval that matches the modelName and savestring
% 4. ./Parameters - lowest fval that matches the modelName
% 5. ./Parameters - lowest fval

if nargin >=3
    saveString = varargin{1};
end
ResultLocation=['.', filesep, 'Results', filesep, resultFolderName];

if ~isempty(dir(sprintf('%s%s*.mat', ResultLocation,filesep))) %if location is not empty, find the best parameter set in the location folder
    saveFolder = sprintf('%s', ResultLocation);
    searchPattern = []; % an emtpy search pattern mean get the best results from the resultLocation folder.
elseif ~isempty(dir(['.', filesep, 'Parameters', filesep, resultFolderName])) % if the results folder is empty check if a coresponding resutls folder exist in the parameters Folder
    saveFolder = ['.', filesep, 'Parameters', filesep, resultFolderName];
    searchPattern = [];
else % if the Results folder is empty use the best parameter vector from the Parameters Folder
    saveFolder = 'Parameters';
    searchPattern = sprintf('%s_FitTo_%s',modelName, strjoin(saveString,'_')); % if the search is limited to the bse parameters folder add a search pattern
    for i = 1:2
        if isempty(findBestSavedParameters(saveFolder,searchPattern)) % if the search pattern yeilds no results change the search pattern
            switch i
                case 1 % if the search pattern with the full save string gives no hits use the modelName as the next search pattern
                    searchPattern = modelName;
                case 2 % if the modelName gives no hits use an empty search pattern i.e. get the best results from the base parameters folder.
                    searchPattern = [];
            end
        else
            break
        end
    end
end

bestParamFileName = findBestSavedParameters(saveFolder,searchPattern); % finde the save file that corresponds to the best solution in the folder saveFolder
load(sprintf('%s%s%s', saveFolder,filesep,bestParamFileName),'theta') % load the veariable theta from the best Parameter Save file
end
