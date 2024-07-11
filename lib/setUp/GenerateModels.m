function []=GenerateModels()
addpath(genpath(fullfile(pwd,'./Models')))

%% find where the AMICI toolbox is  located and add it to path
nStepsOut = 1;% set the number of folders we have to acend to find the "Toolboxes" folder
while and(nStepsOut <= nnz(pwd==filesep),~exist([repmat('../',1,nStepsOut),'Toolboxes/AMICI-0.10.11_SS_eventFix/matlab/'],'dir'))% if the "toolboxes foler doe not exist nStepsOut from the starting folder and we are not at the root folder i.e. we have not earched all folders on the current path.
nStepsOut = nStepsOut+1; % take oun more step out
end
if nStepsOut > nnz(pwd==filesep) %if the number of folder searched excceds the numer of folders in the path return an error'
    error('AMICI toolbox not found! please ensure that the directory %s is avilable from a parent direcory to the current directory','Toolboxes/AMICI-0.10.11_SS_eventFix/' )
end

run([repmat('../',1,nStepsOut),'Toolboxes/AMICI-0.10.11_SS_eventFix/matlab/installAMICI.m']) % set path for the AMICI toolbox.

basePath =pwd;

%%
preExistingModels = dir([repmat('../',1,nStepsOut),'Toolboxes/AMICI-0.10.11_SS_eventFix/models']); %find previously compiled models 
preExistingModels  = preExistingModels(3:end);
preExistingModelNames = vertcat({preExistingModels.name});

symsFiles = vertcat(dir('./Models/*/*_syms.m'),dir('*_syms.m')); % find which _syms.m files exists in all subfolders and current folder.
symsFileNames = vertcat({symsFiles.name})';
modelNames = cellfun(@strrep,symsFileNames,repmat({'_syms.m'},size(symsFileNames)),repmat({''},size(symsFileNames)),'uni',0); % get the model name without the _syms.m surfix

%% find which model are already exists as .mex* files 

mexFiles = vertcat(dir('./Models/*/*.mex*'),dir('*.mex*')); % find which mex files exists in current folder and subfolders.
mexFileNames = vertcat({mexFiles.name})';

mexFileExistsIndx = false(1,size(modelNames,1));
mexIndx = cellfun(@strfind,mexFileNames,repmat({'.mex'},size(mexFileNames)))-1; % find the index of the .mex* surfix

for i = 1:length(mexFileNames) % find the logical idexing for the models that already have existing .mex* files
mexFileExistsIndx(strcmp(mexFileNames{i}(5:mexIndx(i)),modelNames))=true;
end

%% Compile the models that don't have an existing .mex* file  

for i = find(~mexFileExistsIndx)
    if any(strcmp(modelNames{i},preExistingModelNames)) % check if the model name maches any ofh the previousl compiled models 
        modelIdx = strcmp(modelNames{i},preExistingModelNames); % find the index of the previously compiled model that maches the new model 
        rmdir([preExistingModels(modelIdx).folder,filesep,preExistingModels(modelIdx).name],'s') % remove the folder for the previous model as to not interfere with the new model 
    end
   fprintf('\n\nCompiling model: %s\n',modelNames{i});

   if ~exist('./Models/mex_and_simFiles/','dir') % check if a "mex_and_simFiles" directory exists 
       mkdir('./Models/mex_and_simFiles/') % if it does not exist create it.
   end
   try
       amiwrap(modelNames{i},symsFileNames{i}(1:end-2),'./Models/mex_and_simFiles/') % compile the models into the "mex_and_simFiles" folder.
       fprintf('\nModel %s sucessfully compiled.\n',modelNames{i});
   catch ErrorReport
       disp(getReport(ErrorReport));
       fprintf('\nModel %s compilation failed!.\n',modelNames{i});
   end
end

cd(basePath);
end 


