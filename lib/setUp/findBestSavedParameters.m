function [fileName] = findBestSavedParameters(Folder,searchPattern)
% takes a string containitns a folder path and returns a string with the
% file name of the file wilth the lowest number in the file name.
allFiles = dir([Folder, filesep, '*.mat']);
if isempty(allFiles)
    fileName=[];
    return
elseif numel(allFiles)==1
    fileName = allFiles.name;
    return
end 
allFileNames ={allFiles.name};

%%
if nargin == 2 && ~isempty(searchPattern) % if two input arguments are provided and the argument seatch partern is not empty
    if any(contains(allFileNames,searchPattern)) % if any of the files contain the search pattern
        allFileNames = allFileNames(contains(allFileNames,searchPattern)); %consider only the files that match the search pattern
    else %if no files match the search pattern return an empty file name 
        fileName=[];
        return
    end
end

%%
parenthesesIdx = regexp(allFileNames,'[()]');
fValues = cellfun(@(s,i)str2double(s(i(1)+1:i(2)-1)),allFileNames,parenthesesIdx);% extract all cost values in a single vector
[~,minIdx] = min(fValues);%find the idex of the minimum fvalue
fileName = allFileNames{minIdx};

end

