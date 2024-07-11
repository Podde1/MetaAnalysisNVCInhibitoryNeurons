function [] = ExtractMCMC(fileName)

ds = datastore(fileName);
allData = tall(ds);
costData = table2array(allData(:,1));
paramData = table2array(allData(:,2:end-1)); % remove empty column

dTime = regexp(fileName,'\d{6}_\d{6}','match','once');

modelName = regexp(fileName,'Sten2023_v\d{2}','match','once');
if contains(modelName, '_v2')
    ExperimentIdx = [1 4 5 8 9:15]; 
else
    ExperimentIdx = 22:27;
end
variables2Fit = {'HbO','HbT','HbR','CBF','BOLD'};

% min/max cost
fprintf('--- Finding best params ---\n')
MCMCRes = [];
[val, idx] = min(gather(costData));
MCMCRes.minCost.cost = val;
MCMCRes.minCost.params = gather(paramData(idx,:));

[val, idx] = max(gather(costData));
MCMCRes.maxCost.cost = val;
MCMCRes.maxCost.params = gather(paramData(idx,:));

theta = MCMCRes.minCost.params;

% get p names
[~, ~, model, ~, ~, ~, parameterMappingTable, ~] = modelSetup(modelName,ExperimentIdx,variables2Fit);

% configure the total parameter estimation vector to contian different parameter sets for different experiments 
for i = 1:size(parameterMappingTable,1) % loop throught the parameter Mapping table.
    if any(parameterMappingTable.TotalParameterIndex{i}>numel(theta)) % if any of the TotalParameterIndex on th i:th row exceed the number of parameters in the loaded theta-vector. 
        model.parameters(parameterMappingTable.TotalParameterIndex{i}) = model.parameters(parameterMappingTable.ModelParameterIndex{i}); % set the corresponding parameters names for the added parameters.
    end

    if any(parameterMappingTable.TotalParameterIndex{i}>numel(model.parameters)) % if any of the TotalParameterIndex on th i:th row exceed the number of parameters to optimize 
        model.parameters(parameterMappingTable.TotalParameterIndex{i}) = model.parameters(parameterMappingTable.ModelParameterIndex{i}); % set the corresponding parameters names for the added parameters.
    end
end

%% max min parametervalues
fprintf('--- Finding max/min parameter values ---\n')
numParams = gather(size(paramData,2));
for i = 1:numParams
    fprintf('--- Parameter %d out of %d --- \n', i, numParams)
    paramColumn = gather(paramData(:,i));
    [val, idx] = min(paramColumn);
    MCMCRes.(model.parameters{i}).minVal = val;
    MCMCRes.(model.parameters{i}).minValParams = gather(paramData(idx,:));

    [val, idx] = max(paramColumn);
    MCMCRes.(model.parameters{i}).maxVal = val;
    MCMCRes.(model.parameters{i}).maxValParams = gather(paramData(idx,:));
end

%save res
save(sprintf('Results/MCMCRes_%s_%s.mat', modelName, dTime), 'MCMCRes')

filename = sprintf('Results/fval_(%.2f)_%s_%s.mat', MCMCRes.minCost.cost, datestr(now,'yymmdd_HHMMSS'), modelName);
save(filename, 'theta')


%% plotting 
nRows = floor(sqrt(numParams));
nCols =  ceil(sqrt(numParams));
if nRows*nCols < numParams
    nCols = nCols + 1;
end

t = tiledlayout(nRows,nCols);
t.TileSpacing = 'compact';

for i = 1:numParams
    nexttile
    histogram(gather(paramData(:,i)))
    title(model.parameters{i}, 'Interpreter','none')
end

end