function [compatibility] = checkCompatibility()
%% print requirements 
requirementStrings = {'Ensure that a MATLAB release of R2023a or newer is used to run this script.';
'Ensure that the AMICI toolbox is installed  before runnig this script (https://github.com/AMICI-dev/AMICI).';
'Ensure that all relevant models have been compiled before running this script.  To do this please run the "GenerateModels()" script in MATLAB R2017b.'};

nRequirements = numel(requirementStrings);

fprintf('Please ensure that the following %d requirements are met before plotting the results: \n', nRequirements)

for i = 1:nRequirements
    fprintf('%d. %s \n',i,requirementStrings{i});
end
%% check requirments
% check matlab verision 
versionStr = version;
versionIdx = regexp(version,'\(R20\d\d[a-b]\)');

% check existance of mex files 
mexfiles=dir('.\Models\mex_and_simFiles\*.mex*');


logicalChecks = [str2double(versionStr(versionIdx+[4:5])) < 23;     % check matlab verision 
                 ~contains(path,'AMICI');                           % check if "AMICI" is in path
                 isempty(mexfiles)];                                % check existance of mex files 

warningMsg ={sprintf('Plotting scripts require a MATLAB release of R2023a or newer. Current version is: %s! ',version);
            'AMICI is not found in matlab''s search path!'
            'No ".mex"-files were found in the ".\Models\mex_and_simFiles\"-directory.'};

if any(logicalChecks)
    warningIdx = [1:nRequirements];
    for i = warningIdx(logicalChecks)
        warning('\n%s Requirement %d. might not be fulfilled!', warningMsg{i},i)
    end
end

%% have used verify requirments 
fprintf('\n \n');

check = 'n';

while ~strcmp(check,'y')
    check = input('Are all the above requirements fulfilled? [y/n] >> ','s');

    fprintf('\n');

    if and(~strcmp(check,'y'),~strcmp(check,'n'))
        fprintf('\n Your answer must be y - YES or n - NO. %s is not a valid answer!\n', check)
    end
    if strcmp(check,'n')
        fprintf('User deems that the requirements NOT fulfilled! Terminating script.\n')
        compatibility = 0;
        return
    end
end

fprintf('User verifies that the requirements are fulfilled! Plotting reults. \n')
compatibility = 1;
end 