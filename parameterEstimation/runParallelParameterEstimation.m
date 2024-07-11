function [] = runParallelParameterEstimation(resultFolderName, modelName, ExperimentIdx, variables2Fit, nTimePerRun, totalMaxRunTime)

fprintf('Estimated to finnish by %s\n',datestr(now+hours(fix(totalMaxRunTime/nTimePerRun)*nTimePerRun) + seconds(90))) % print estimated finish time based on totalRunTime and an estimated time it takes to start the local parallel pool 

parPoolObj = parpool('local'); % start local parallel pool
nWorkers = parPoolObj.NumWorkers; % get the numeber of workers in the parallel pool 

parfor i = 1:nWorkers*(fix(totalMaxRunTime/nTimePerRun)) % pun Parameter estimation equal to the numebr of workers*numer of whole runs TimePerRun fits into total run time 
    parameterEstimation(resultFolderName, modelName, ExperimentIdx, variables2Fit, nTimePerRun);
end

delete(parPoolObj); % close local parallel pool