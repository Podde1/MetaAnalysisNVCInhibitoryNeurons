function [model] = getModelNames(modelsymFileName)
%% takes the file name of an AMICI formulated model amnd returns teh snames of the model states, parameters, constatns and observables

    rawStr = readcell([modelsymFileName,'.m'],'FileType','text','Delimiter',';');  %read the model file as a text file
    rowIdentifierStr = {'model.sym.x = ','model.sym.p = ','model.sym.k = ','model.sym.y = '}; % pre define the identifyerstring used to find the rows with the names
    modelfieldNames={'states', 'parameters', 'constants', 'observables'}; % pre define the output structure fieldnames
    x = cell(1,4); % pre allocate a cell array to store the partally formated strngs

    for i = 1:4 % loop over four itterations 
        x{i} = strsplit(rawStr{startsWith(rawStr(:,1),rowIdentifierStr{i}),1},{'[',',',' ',']'});% find the row that starts with "model.sym.X = ". and split the string at various delimiters. 
        x{i} = x{i}(find(strcmp(x{i},'='))+1:end); % remove the "model.sym.x" and "=" cells from the list.

        % assign the strings to the model fields. The ~strcmp command removes any empty cells
        model.(modelfieldNames{i}) = x{i}(~strcmp(x{i},''));
    end
end